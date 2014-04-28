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

class IGIIntegrator extends SurfaceIntegrator {
  IGIIntegrator(int nl, int ns, double rrt, int maxd, double gl, int ng) {
    nLightPaths = RoundUpPow2(nl);
    nLightSets = RoundUpPow2(ns);
    rrThreshold = rrt;
    maxSpecularDepth = maxd;
    virtualLights = new List(nLightSets);
    gLimit = gl;
    nGatherSamples = ng;

    for (int i = 0; i < nLightSets; ++i) {
      virtualLights[i] = new List<_VirtualLight>();
    }
  }

  Spectrum Li(Scene scene, Renderer renderer,
      RayDifferential ray, Intersection isect,
      Sample sample, RNG rng) {
    Spectrum L = new Spectrum(0.0);
    Vector wo = -ray.direction;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Evaluate BSDF at hit point
    BSDF bsdf = isect.getBSDF(ray);
    Point p = bsdf.dgShading.p;
    Normal n = bsdf.dgShading.nn;

    L += Integrator.UniformSampleAllLights(scene, renderer, p, n,
                                wo, isect.rayEpsilon, ray.time, bsdf, sample,
                                rng, lightSampleOffsets, bsdfSampleOffsets);

    // Compute indirect illumination with virtual lights
    int lSet = Math.min((sample.oneD[vlSetOffset][0] * nLightSets).toInt(),
                        nLightSets - 1);

    for (int i = 0; i < virtualLights[lSet].length; ++i) {
      _VirtualLight vl = virtualLights[lSet][i];

      // Compute virtual light's tentative contribution _Llight_
      double d2 = Vector.DistanceSquared(p, vl.p);
      Vector wi = Vector.Normalize(vl.p - p);
      double G = Vector.AbsDot(wi, n) * Vector.AbsDot(wi, vl.n) / d2;

      G = Math.min(G, gLimit);

      Spectrum f = bsdf.f(wo, wi);

      if (G == 0.0 || f.isBlack()) {
        continue;
      }

      Spectrum Llight = f * G * vl.pathContrib / nLightPaths;

      RayDifferential connectRay = new RayDifferential.child(p, wi, ray,
                                        isect.rayEpsilon,
                                        Math.sqrt(d2) * (1.0 - vl.rayEpsilon));

      Llight *= renderer.transmittance(scene, connectRay, null, rng);

      // Possibly skip virtual light shadow ray with Russian roulette
      if (Llight.luminance() < rrThreshold) {
        double continueProbability = 0.1;
        if (rng.randomFloat() > continueProbability) {
          continue;
        }

        Llight /= continueProbability;
      }

      // Add contribution from _VirtualLight_ _vl_
      if (!scene.intersectP(connectRay)) {
        L += Llight;
      }
    }

    if (ray.depth < maxSpecularDepth) {
      // Do bias compensation for bounding geometry term
      int nSamples = (ray.depth == 0) ? nGatherSamples : 1;
      List<double> pdf = [0.0];

      for (int i = 0; i < nSamples; ++i) {
        Vector wi = new Vector();
        BSDFSample bsdfSample = (ray.depth == 0) ?
                new BSDFSample.sample(sample, gatherSampleOffset, i) :
                new BSDFSample.random(rng);
                Spectrum f = bsdf.sample_f(wo, wi, bsdfSample, pdf,
                                   BSDF_ALL & ~BSDF_SPECULAR);

        if (!f.isBlack() && pdf[0] > 0.0) {
          // Trace ray for bias compensation gather sample
          double maxDist = Math.sqrt(Vector.AbsDot(wi, n) / gLimit);
          RayDifferential gatherRay = new RayDifferential.child(p, wi, ray,
                                                    isect.rayEpsilon, maxDist);

          Intersection gatherIsect = new Intersection();
          Spectrum Li = renderer.Li(scene, gatherRay, sample, rng, gatherIsect);
          if (Li.isBlack()) {
            continue;
          }

          // Add bias compensation ray contribution to radiance sum
          double Ggather = Vector.AbsDot(wi, n) *
                           Vector.AbsDot(-wi, gatherIsect.dg.nn) /
                           Vector.DistanceSquared(p, gatherIsect.dg.p);
          if (Ggather - gLimit > 0.0 && Ggather.isFinite) {
            double gs = (Ggather - gLimit) / Ggather;
            L += f * Li * (Vector.AbsDot(wi, n) * gs / (nSamples * pdf[0]));
          }
        }
      }
    }

    if (ray.depth + 1 < maxSpecularDepth) {
      Vector wi = new Vector();
      // Trace rays for specular reflection and refraction
      L += Integrator.SpecularReflect(ray, bsdf, rng, isect, renderer, scene,
                                      sample);
      L += Integrator.SpecularTransmit(ray, bsdf, rng, isect, renderer, scene,
                                       sample);
    }

    return L;
  }

  void requestSamples(Sampler sampler, Sample sample, Scene scene) {
    // Allocate and request samples for sampling all lights
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

    vlSetOffset = sample.add1D(1);

    if (sampler != null) {
      nGatherSamples = sampler.roundSize(nGatherSamples);
    }

    gatherSampleOffset = new BSDFSampleOffsets(nGatherSamples, sample);
  }

  void preprocess(Scene scene, Camera camera, Renderer renderer) {
    if (scene.lights.isEmpty) {
      return;
    }

    RNG rng = new RNG();

    // Compute samples for emitted rays from lights
    Float32List lightNum = new Float32List(nLightPaths * nLightSets);
    Float32List lightSampPos = new Float32List(2 * nLightPaths * nLightSets);
    Float32List lightSampComp = new Float32List(nLightPaths * nLightSets);
    Float32List lightSampDir = new Float32List(2 * nLightPaths * nLightSets);

    LDShuffleScrambled1D(nLightPaths, nLightSets, lightNum, rng);
    LDShuffleScrambled2D(nLightPaths, nLightSets, lightSampPos, rng);
    LDShuffleScrambled1D(nLightPaths, nLightSets, lightSampComp, rng);
    LDShuffleScrambled2D(nLightPaths, nLightSets, lightSampDir, rng);

    // Precompute information for light sampling densities
    Distribution1D lightDistribution = Integrator.ComputeLightSamplingCDF(scene);

    List<double> lightPdf = [0.0];
    List<double> pdf = [0.0];

    for (int s = 0; s < nLightSets; ++s) {
      for (int i = 0; i < nLightPaths; ++i) {
        // Follow path _i_ from light to create virtual lights
        int sampOffset = s * nLightPaths + i;

        // Choose light source to trace virtual light path from
        int ln = lightDistribution.sampleDiscrete(lightNum[sampOffset],
                                                  lightPdf);
        Light light = scene.lights[ln];

        // Sample ray leaving light source for virtual light path
        RayDifferential ray = new RayDifferential();

        LightSample ls = new LightSample(lightSampPos[2 * sampOffset],
                                         lightSampPos[2 * sampOffset + 1],
                                         lightSampComp[sampOffset]);

        Normal Nl = new Normal();
        Spectrum alpha = light.sampleL(scene, ls, lightSampDir[2 * sampOffset],
                                       lightSampDir[2 * sampOffset + 1],
                                       camera.shutterOpen, ray, Nl, pdf);

        if (pdf[0] == 0.0 || alpha.isBlack()) {
          continue;
        }

        alpha = alpha * Vector.AbsDot(Nl, ray.direction) /
                (pdf[0] * lightPdf[0]);

        Intersection isect = new Intersection();
        while (scene.intersect(ray, isect) && !alpha.isBlack()) {
          // Create virtual light and sample new ray for path
          alpha *= renderer.transmittance(scene,
                                          new RayDifferential.fromRay(ray),
                                          null, rng);

          Vector wo = -ray.direction;
          BSDF bsdf = isect.getBSDF(ray);

          // Create virtual light at ray intersection point
          Spectrum contrib = alpha * bsdf.rho2(wo, rng) / Math.PI;

          virtualLights[s].add(new _VirtualLight(isect.dg.p, isect.dg.nn,
                                                 contrib, isect.rayEpsilon));

          // Sample new ray direction and update weight for virtual light path
          Vector wi = new Vector();
          BSDFSample bsdfSample = new BSDFSample.random(rng);
          Spectrum fr = bsdf.sample_f(wo, wi, bsdfSample, pdf);
          if (fr.isBlack() || pdf == 0.0) {
            break;
          }

          Spectrum contribScale = fr * Vector.AbsDot(wi, bsdf.dgShading.nn) /
                                  pdf[0];

          // Possibly terminate virtual light path with Russian roulette
          double rrProb = Math.min(1.0, contribScale.luminance());
          if (rng.randomFloat() > rrProb) {
            break;
          }

          alpha *= contribScale / rrProb;

          ray = new RayDifferential.child(isect.dg.p, wi, ray,
                                               isect.rayEpsilon);
        }
      }
    }
  }

  static IGIIntegrator Create(ParamSet params) {
    int nLightPaths = params.findOneInt("nlights", 64);
    int nLightSets = params.findOneInt("nsets", 4);
    double rrThresh = params.findOneFloat("rrthreshold", 0.0001);
    int maxDepth = params.findOneInt("maxdepth", 5);
    double glimit = params.findOneFloat("glimit", 10.0);
    int gatherSamples = params.findOneInt("gathersamples", 16);

    return new IGIIntegrator(nLightPaths, nLightSets, rrThresh,
                             maxDepth, glimit, gatherSamples);
  }

  List<LightSampleOffsets> lightSampleOffsets;
  List<BSDFSampleOffsets> bsdfSampleOffsets;
  int nLightPaths, nLightSets;
  double gLimit;
  int nGatherSamples;
  double rrThreshold;
  int maxSpecularDepth;
  int vlSetOffset;
  BSDFSampleOffsets gatherSampleOffset;
  List<List<_VirtualLight>> virtualLights;
}

class _VirtualLight {
  _VirtualLight(this.p, this.n, this.pathContrib, this.rayEpsilon);

  Point p;
  Normal n;
  Spectrum pathContrib;
  double rayEpsilon;
}
