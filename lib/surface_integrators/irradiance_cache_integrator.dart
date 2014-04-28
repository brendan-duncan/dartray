/****************************************************************************
 * Copyright (C) 2014 by Brendan Duncan.                                    *
 *                                                                          *
 * This file is part of DartRay.                                            *
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License");          *
 * you may not use this file except in compliance with the License.         *
 * You may obtain a copy of the License at                                  *
 *                                                                          *
 * http://www.apache.org/licenses/LICENSE-2.0                               *
 *                                                                          *
 * Unless required by applicable law or agreed to in writing, software      *
 * distributed under the License is distributed on an "AS IS" BASIS,        *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 * See the License for the specific language governing permissions and      *
 * limitations under the License.                                           *
 *                                                                          *
 * This project is based on PBRT v2 ; see http://www.pbrt.org               *
 * pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.  *
 ****************************************************************************/
part of surface_integrators;

class IrradianceCacheIntegrator extends SurfaceIntegrator {
  IrradianceCacheIntegrator(this.minWeight, this.minSamplePixelSpacing,
                            this.maxSamplePixelSpacing,
                            double maxang, this.maxSpecularDepth,
                            this.maxIndirectDepth, this.nSamples) {
    cosMaxSampleAngleDifference = Math.cos(Degrees(maxang));
  }

  Spectrum Li(Scene scene, Renderer renderer, RayDifferential ray,
              Intersection isect, Sample sample, RNG rng) {
    Spectrum L = new Spectrum(0.0);

    // Evaluate BSDF at hit point
    BSDF bsdf = isect.getBSDF(ray);
    Vector wo = -ray.direction;
    Point p = bsdf.dgShading.p;
    Normal n = bsdf.dgShading.nn;
    L += isect.Le(wo);

    // Compute direct lighting for irradiance cache
    L += Integrator.UniformSampleAllLights(scene, renderer, p, n, wo,
             isect.rayEpsilon, ray.time, bsdf, sample, rng,
             lightSampleOffsets, bsdfSampleOffsets);

    // Compute indirect lighting for irradiance cache
    if (ray.depth + 1 < maxSpecularDepth) {
      // Trace rays for specular reflection and refraction
      L += Integrator.SpecularReflect(ray, bsdf, rng, isect, renderer, scene,
                                      sample);

      L += Integrator.SpecularTransmit(ray, bsdf, rng, isect, renderer, scene,
                                       sample);
    }

    // Estimate indirect lighting with irradiance cache
    Normal ng = Normal.FaceForward(isect.dg.nn, wo);

    // Compute pixel spacing in world space at intersection point
    double pixelSpacing = Math.sqrt(Vector.Cross(isect.dg.dpdx,
                                                 isect.dg.dpdy).length());

    int flags = BSDF_REFLECTION | BSDF_DIFFUSE | BSDF_GLOSSY;
    L += indirectLo(p, ng, pixelSpacing, wo, isect.rayEpsilon,
                    bsdf, flags, rng, scene, renderer);

    flags = BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY;
    L += indirectLo(p, -ng, pixelSpacing, wo, isect.rayEpsilon,
                    bsdf, flags, rng, scene, renderer);

    return L;
  }

  void requestSamples(Sampler sampler, Sample sample, Scene scene) {
    int nLights = scene.lights.length;
    lightSampleOffsets = new List<LightSampleOffsets>(nLights);
    bsdfSampleOffsets = new List<BSDFSampleOffsets>(nLights);
    for (int i = 0; i < nLights; ++i) {
      Light light = scene.lights[i];
      int nSamples = light.nSamples;
      if (sampler != null) {
        nSamples = sampler.roundSize(nSamples);
      }
      lightSampleOffsets[i] = new LightSampleOffsets(nSamples, sample);
      bsdfSampleOffsets[i] = new BSDFSampleOffsets(nSamples, sample);
    }
  }

  void preprocess(Scene scene, Camera camera, Renderer renderer) {
    BBox wb = new BBox.from(scene.worldBound);
    Vector delta = (wb.pMax - wb.pMin) * 0.01;
    wb.pMin -= delta;
    wb.pMax += delta;

    octree = new Octree(wb);

    // Prime irradiance cache
    minWeight *= 1.5;

    List<int> extent = [0, 0, 0, 0];
    camera.film.getSampleExtent(extent);

    HaltonSampler sampler = new HaltonSampler(extent[0], extent[1],
                                              extent[2], extent[3],
                                              camera.shutterOpen,
                                              camera.shutterClose, 1);

    Sample sample = new Sample(sampler, this, null, scene);

    IrradiancePrimeTask task = new IrradiancePrimeTask(scene, renderer, camera,
                                                       sampler, sample, this);
    task.run();

    minWeight /= 1.5;
  }

  Spectrum indirectLo(Point p, Normal ng, double pixelSpacing,
                      Vector wo, double rayEpsilon, BSDF bsdf,
                      int flags, RNG rng, Scene scene, Renderer renderer) {
    if (bsdf.numComponents(flags) == 0) {
      return new Spectrum(0.0);
    }

    Spectrum E = new Spectrum();
    Vector wi = new Vector();

    // Get irradiance _E_ and average incident direction _wi_ at point _p_
    if (!interpolateE(scene, p, ng, E, wi)) {
      // Compute irradiance at current point
      Stats.IRRADIANCE_CACHE_STARTED_COMPUTING_IRRADIANCE(p, ng);
      List<int> scramble = [rng.randomUint(), rng.randomUint()];
      double minHitDistance = INFINITY;
      Vector wAvg = new Vector(0.0, 0.0, 0.0);
      Spectrum LiSum = new Spectrum(0.0);

      List<double> u = [0.0, 0.0];
      for (int i = 0; i < nSamples; ++i) {
        // Sample direction for irradiance estimate ray
        Sample02(i, scramble, u);

        Vector w = CosineSampleHemisphere(u[0], u[1]);

        RayDifferential r = new RayDifferential(p, bsdf.localToWorld(w),
                                                rayEpsilon);
        r.direction = Vector.FaceForward(r.direction, ng);

        // Trace ray to sample radiance for irradiance estimate
        Stats.IRRADIANCE_CACHE_STARTED_RAY(r);
        Spectrum L = pathL(r, scene, renderer, rng);
        LiSum += L;
        wAvg += r.direction * L.luminance();
        minHitDistance = Math.min(minHitDistance, r.maxDistance);
        Stats.IRRADIANCE_CACHE_FINISHED_RAY(r, r.maxDistance, L);
      }

      E = LiSum * (Math.PI / nSamples);
      Stats.IRRADIANCE_CACHE_FINISHED_COMPUTING_IRRADIANCE(p, ng);

      // Add computed irradiance value to cache

      // Compute irradiance sample's contribution extent and bounding box
      double maxDist = maxSamplePixelSpacing * pixelSpacing;
      double minDist = minSamplePixelSpacing * pixelSpacing;
      double contribExtent = (minHitDistance / 2.0).clamp(minDist, maxDist);

      BBox sampleExtent = new BBox(p);
      sampleExtent.expand(contribExtent);
      Stats.IRRADIANCE_CACHE_ADDED_NEW_SAMPLE(p, ng, contribExtent, E, wAvg,
                                              pixelSpacing);

      // Allocate _IrradianceSample_, get write lock, add to octree
      IrradianceSample sample = new IrradianceSample(E, p, ng, wAvg,
                                                     contribExtent);
      octree.add(sample, sampleExtent);

      wi = wAvg;
    }

    // Compute reflected radiance due to irradiance and BSDF
    if (wi.lengthSquared() == 0.0) {
      return new Spectrum(0.0);
    }

    return bsdf.f(wo, Vector.Normalize(wi), flags) * E;
  }

  bool interpolateE(Scene scene, Point p, Normal n, Spectrum E, Vector wi) {
    if (octree == null) {
      return false;
    }

    Stats.IRRADIANCE_CACHE_STARTED_INTERPOLATION(p, n);
    IrradProcess proc = new IrradProcess(p, n, minWeight,
                                         cosMaxSampleAngleDifference);

    octree.lookup(p, proc);

    Stats.IRRADIANCE_CACHE_FINISHED_INTERPOLATION(p, n,
                                                  proc.successful() ? 1 : 0,
                                                  proc.nFound);

    if (!proc.successful()) {
      return false;
    }

    E.copy(proc.getIrradiance());
    wi.copy(proc.getAverageDirection());

    return true;
  }

  Spectrum pathL(Ray r, Scene scene, Renderer renderer, RNG rng) {
    Spectrum L = new Spectrum(0.0);
    Spectrum pathThroughput = new Spectrum(1.0);
    RayDifferential ray = new RayDifferential.fromRay(r);
    bool specularBounce = false;

    for (int pathLength = 0; ; ++pathLength) {
      // Find next vertex of path
      Intersection isect = new Intersection();
      if (!scene.intersect(ray, isect)) {
        break;
      }

      if (pathLength == 0) {
        r.maxDistance = ray.maxDistance;
      }

      pathThroughput *= renderer.transmittance(scene, ray, null, rng);

      // Possibly add emitted light at path vertex
      if (specularBounce) {
        L += pathThroughput * isect.Le(-ray.direction);
      }

      // Evaluate BSDF at hit point
      BSDF bsdf = isect.getBSDF(ray);

      // Sample illumination from lights to find path contribution
      Point p = bsdf.dgShading.p;
      Normal n = bsdf.dgShading.nn;
      Vector wo = -ray.direction;

      L += pathThroughput *
            Integrator.UniformSampleOneLight(scene, renderer, p, n, wo,
                                             isect.rayEpsilon, ray.time, bsdf,
                                             null, rng);

      if (pathLength + 1 == maxIndirectDepth) {
        break;
      }

      // Sample BSDF to get new path direction
      // Get random numbers for sampling new direction, bs1, bs2, and bcs.
      Vector wi = new Vector();
      List<double> pdf = [0.0];
      List<int> flags = [0];

      Spectrum f = bsdf.sample_f(wo, wi, new BSDFSample.random(rng),
                                 pdf, BSDF_ALL, flags);

      if (f.isBlack() || pdf[0] == 0.0) {
        break;
      }

      specularBounce = (flags[0] & BSDF_SPECULAR) != 0;
      pathThroughput *= f * Vector.AbsDot(wi, n) / pdf[0];

      ray = new RayDifferential.child(p, wi, ray, isect.rayEpsilon);

      // Possibly terminate the path
      if (pathLength > 2) {
        double rrProb = Math.min(1.0, pathThroughput.luminance());

        if (rng.randomFloat() > rrProb) {
          break;
        }

        pathThroughput /= rrProb;
      }
    }

    return L;
  }

  static IrradianceCacheIntegrator Create(ParamSet params) {
    double minWeight = params.findOneFloat('minweight', 0.5);
    double minSpacing = params.findOneFloat('minpixelspacing', 2.5);
    double maxSpacing = params.findOneFloat('maxpixelspacing', 15.0);
    double maxAngle = params.findOneFloat('maxangledifference', 10.0);
    int maxSpecularDepth = params.findOneInt('maxspeculardepth', 5);
    int maxIndirectDepth = params.findOneInt('maxindirectdepth', 3);
    int nSamples = params.findOneInt('nsamples', 4096);

    return new IrradianceCacheIntegrator(minWeight, minSpacing, maxSpacing,
                                         maxAngle, maxSpecularDepth,
                                         maxIndirectDepth, nSamples);
  }

  double minSamplePixelSpacing;
  double maxSamplePixelSpacing;
  double minWeight;
  double cosMaxSampleAngleDifference;
  int nSamples;
  int maxSpecularDepth;
  int maxIndirectDepth;

  // Declare sample parameters for light source sampling
  List<LightSampleOffsets> lightSampleOffsets;
  List<BSDFSampleOffsets> bsdfSampleOffsets;
  Octree octree;
}

class IrradianceSample {
  IrradianceSample(this.E, this.p, this.n, this.wAvg, this.maxDist);

  Spectrum E;
  Normal n;
  Point p;
  Vector wAvg;
  double maxDist;
}

class IrradiancePrimeTask {
  IrradiancePrimeTask(this.scene, this.renderer, this.camera, this.sampler,
                      this.origSample, this.irradianceCache);

  void run() {
    if (sampler == null) {
      return;
    }

    int sampleCount;
    RNG rng = new RNG(29);
    int maxSamples = sampler.maximumSampleCount();
    List<Sample> samples = origSample.duplicate(maxSamples);

    while ((sampleCount = sampler.getMoreSamples(samples, rng)) > 0) {
      for (int i = 0; i < sampleCount; ++i) {
        RayDifferential ray = new RayDifferential();
        camera.generateRayDifferential(samples[i], ray);

        Intersection isect = new Intersection();
        if (scene.intersect(ray, isect)) {
          irradianceCache.Li(scene, renderer, ray, isect, samples[i], rng);
        }
      }
    }
  }

  Scene scene;
  Camera camera;
  Renderer renderer;
  Sampler sampler;
  Sample origSample;
  IrradianceCacheIntegrator irradianceCache;
}

class IrradProcess {
  IrradProcess(this.p, this.n, this.minWeight,
               this.cosMaxSampleAngleDifference) {
    nFound = 0;
    sumWt = 0.0;
    E = new Spectrum(0.0);
    wAvg = new Vector(0.0, 0.0, 0.0);
  }

  bool call(IrradianceSample sample) {
    // Compute estimate error term and possibly use sample
    double perr = Vector.Distance(p, sample.p) / sample.maxDist;
    double nerr = Math.sqrt((1.0 - Vector.Dot(n, sample.n)) /
                            (1.0 - cosMaxSampleAngleDifference));
    double err = Math.max(perr, nerr);
    Stats.IRRADIANCE_CACHE_CHECKED_SAMPLE(sample, perr, nerr);
    if (err < 1.0) {
      ++nFound;
      double wt = 1.0 - err;
      E += sample.E * wt;
      wAvg += sample.wAvg * wt;
      sumWt += wt;
    }
    return true;
  }

  bool successful() {
    return sumWt >= minWeight;
  }

  Spectrum getIrradiance() {
    return E / sumWt;
  }

  Vector getAverageDirection() {
    return wAvg;
  }

  Point p;
  Normal n;
  double minWeight;
  double cosMaxSampleAngleDifference;
  double sumWt;
  int nFound;
  Spectrum E;
  Vector wAvg;
}
