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
part of core;

/**
 * Numerically integrates a light transport equation.
 *
 * Integrators handle the task of simulating the propogation of light in the
 * scene in order to compute how much light arives at image sample positions
 * on the film plane.
 *
 * Defines the common interface for [SurfaceIntegrator] and [VolumeIntegrator].
 */
abstract class Integrator {
  void preprocess(Scene scene, Camera camera, Renderer renderer) {
  }

  void requestSamples(Sampler sampler, Sample sample, Scene scene) {
  }

  static Spectrum UniformSampleAllLights(Scene scene, Renderer renderer,
                                    Point p, Normal n, Vector wo,
                                    double rayEpsilon, double time,
                                    BSDF bsdf, Sample sample,
                                    RNG rng,
                                    List<LightSampleOffsets> lightSampleOffsets,
                                    List<BSDFSampleOffsets> bsdfSampleOffsets) {
    Spectrum L = new Spectrum(0.0);

    for (int i = 0; i < scene.lights.length; ++i) {
      Light light = scene.lights[i];
      int nSamples = lightSampleOffsets != null ?
                     lightSampleOffsets[i].nSamples : 1;

      // Estimate direct lighting from light samples
      Spectrum Ld = new Spectrum(0.0);

      for (int j = 0; j < nSamples; ++j) {
        // Find light and BSDF sample values for direct lighting estimate
        LightSample lightSample;
        BSDFSample bsdfSample;
        if (lightSampleOffsets != null && bsdfSampleOffsets != null) {
          lightSample = new LightSample.sample(sample, lightSampleOffsets[i], j);
          bsdfSample = new BSDFSample.sample(sample, bsdfSampleOffsets[i], j);
        } else {
          lightSample = new LightSample.random(rng);
          bsdfSample = new BSDFSample.random(rng);
        }

        Ld += EstimateDirect(scene, renderer, light, p, n, wo,
                             rayEpsilon, time, bsdf, rng, lightSample,
                             bsdfSample, BSDF_ALL & ~BSDF_SPECULAR);
      }

      L += Ld / nSamples.toDouble();
    }

    return L;
  }

  static Spectrum UniformSampleOneLight(Scene scene, Renderer renderer,
                                        Point p, Normal n, Vector wo,
                                        double rayEpsilon, double time,
                                        BSDF bsdf, Sample sample, RNG rng,
                                        [int lightNumOffset = -1,
                                        LightSampleOffsets lightSampleOffset,
                                        BSDFSampleOffsets bsdfSampleOffset]) {
    // Randomly choose a single light to sample
    int nLights = scene.lights.length;
    if (nLights == 0) {
      return new Spectrum(0.0);
    }

    int lightNum;
    if (lightNumOffset != -1) {
      lightNum = (sample.oneD[lightNumOffset][0] * nLights).floor();
    } else {
      lightNum = (rng.randomFloat() * nLights).floor();
    }

    lightNum = Math.min(lightNum, nLights - 1);
    Light light = scene.lights[lightNum];

    // Initialize light and bsdf samples for single light sample
    LightSample lightSample;
    BSDFSample bsdfSample;
    if (lightSampleOffset != null && bsdfSampleOffset != null) {
      lightSample = new LightSample.sample(sample, lightSampleOffset, 0);
      bsdfSample = new BSDFSample.sample(sample, bsdfSampleOffset, 0);
    } else {
      lightSample = new LightSample.random(rng);
      bsdfSample = new BSDFSample.random(rng);
    }

    double s = nLights.toDouble();
    return EstimateDirect(scene, renderer, light, p, n, wo,
                          rayEpsilon, time, bsdf, rng, lightSample,
                          bsdfSample, BSDF_ALL & ~BSDF_SPECULAR) * s;
  }

  static Spectrum EstimateDirect(Scene scene, Renderer renderer,
                                 Light light, Point p, Normal n, Vector wo,
                                 double rayEpsilon, double time, BSDF bsdf,
                                 RNG rng, LightSample lightSample,
                                 BSDFSample bsdfSample, int flags) {
    Spectrum Ld = new Spectrum(0.0);
    // Sample light source with multiple importance sampling
    Vector wi = new Vector();
    List<double> lightPdf = [0.0];
    List<double> bsdfPdf = [0.0];
    VisibilityTester visibility = new VisibilityTester();
    Spectrum Li = light.sampleLAtPoint(p, rayEpsilon, lightSample, time,
                                       wi, lightPdf, visibility);

    if (lightPdf[0] > 0.0 && !Li.isBlack()) {
      Spectrum f = bsdf.f(wo, wi, flags);
      if (!f.isBlack() && visibility.unoccluded(scene)) {
        // Add light's contribution to reflected radiance
        Li *= visibility.transmittance(scene, renderer, null, rng);
        if (light.isDeltaLight()) {
          Ld += f * Li * (Vector.AbsDot(wi, n) / lightPdf[0]);
        } else {
          bsdfPdf[0] = bsdf.pdf(wo, wi, flags);
          double weight = PowerHeuristic(1, lightPdf[0], 1, bsdfPdf[0]);
          Ld += f * Li * ((Vector.AbsDot(wi, n) * weight / lightPdf[0]));
        }
      }
    }

    // Sample BSDF with multiple importance sampling
    if (!light.isDeltaLight()) {
      List<int> sampledType = [0];
      Spectrum f = bsdf.sample_f(wo, wi, bsdfSample, bsdfPdf, flags,
                                    sampledType);
      if (!f.isBlack() && bsdfPdf[0] > 0.0) {
        double weight = 1.0;
        if ((sampledType[0] & BSDF_SPECULAR) == 0) {
          lightPdf[0] = light.pdf(p, wi);
          if (lightPdf[0] == 0.0) {
            return Ld;
          }
          weight = PowerHeuristic(1, bsdfPdf[0], 1, lightPdf[0]);
        }

        // Add light contribution from BSDF sampling
        Intersection lightIsect = new Intersection();
        Spectrum Li = new Spectrum(0.0);
        RayDifferential ray = new RayDifferential(p, wi, rayEpsilon, INFINITY,
                                                  time);

        if (scene.intersect(ray, lightIsect)) {
          if (lightIsect.primitive.getAreaLight() == light) {
            Li = lightIsect.Le(-wi);
          }
        } else {
          Li = light.Le(ray);
        }

        if (!Li.isBlack()) {
          Li *= renderer.transmittance(scene, ray, null, rng);
          Ld += f * Li * (Vector.AbsDot(wi, n) * weight / bsdfPdf[0]);
        }
      }
    }

    return Ld;
  }

  static Spectrum SpecularReflect(RayDifferential ray, BSDF bsdf, RNG rng,
                                  Intersection isect, Renderer renderer,
                                  Scene scene, Sample sample) {
    Vector wo = -ray.direction;
    Vector wi = new Vector();
    List<double> pdf = [0.0];
    Point p = bsdf.dgShading.p;
    Normal n = bsdf.dgShading.nn;
    Spectrum f = bsdf.sample_f(wo, wi, new BSDFSample.random(rng), pdf,
                               BSDF_REFLECTION | BSDF_SPECULAR);
    Spectrum L = new Spectrum(0.0);

    if (pdf[0] > 0.0 && !f.isBlack() && Vector.AbsDot(wi, n) != 0.0) {
      // Compute ray differential rd for specular reflection
      RayDifferential rd = new RayDifferential.child(p, wi, ray,
                                                     isect.rayEpsilon);
      if (ray.hasDifferentials) {
        rd.hasDifferentials = true;
        rd.rxOrigin = p + isect.dg.dpdx;
        rd.ryOrigin = p + isect.dg.dpdy;
        // Compute differential reflected directions
        Normal dndx = bsdf.dgShading.dndu * bsdf.dgShading.dudx +
                      bsdf.dgShading.dndv * bsdf.dgShading.dvdx;
        Normal dndy = bsdf.dgShading.dndu * bsdf.dgShading.dudy +
                      bsdf.dgShading.dndv * bsdf.dgShading.dvdy;
        Vector dwodx = -ray.rxDirection - wo;
        Vector dwody = -ray.ryDirection - wo;
        double dDNdx = Vector.Dot(dwodx, n) + Vector.Dot(wo, dndx);
        double dDNdy = Vector.Dot(dwody, n) + Vector.Dot(wo, dndy);

        rd.rxDirection = wi - dwodx + new Vector.from(dndx * Vector.Dot(wo, n) +
                         n * dDNdx) * 2.0;

        rd.ryDirection = wi - dwody + new Vector.from(dndy * Vector.Dot(wo, n) +
                         n * dDNdy)* 2.0;
      }

      Stats.STARTED_SPECULAR_REFLECTION_RAY(rd);
      Spectrum Li = renderer.Li(scene, rd, sample, rng);
      L = f * Li * (Vector.AbsDot(wi, n) / pdf[0]);
      Stats.FINISHED_SPECULAR_REFLECTION_RAY(rd);
    }

    return L;
  }

  static Spectrum SpecularTransmit(RayDifferential ray, BSDF bsdf, RNG rng,
                                   Intersection isect, Renderer renderer,
                                   Scene scene, Sample sample) {
    Vector wo = -ray.direction;
    Vector wi = new Vector();
    List<double> pdf = [0.0];
    Point p = bsdf.dgShading.p;
    Normal n = bsdf.dgShading.nn;
    Spectrum f = bsdf.sample_f(wo, wi, new BSDFSample.random(rng), pdf,
                               BSDF_TRANSMISSION | BSDF_SPECULAR);
    Spectrum L = new Spectrum(0.0);

    if (pdf[0] > 0.0 && !f.isBlack() && Vector.AbsDot(wi, n) != 0.0) {
      // Compute ray differential _rd_ for specular transmission
      RayDifferential rd = new RayDifferential.child(p, wi, ray,
                                                     isect.rayEpsilon);

      if (ray.hasDifferentials) {
        rd.hasDifferentials = true;
        rd.rxOrigin = p + isect.dg.dpdx;
        rd.ryOrigin = p + isect.dg.dpdy;

        double eta = bsdf.eta;
        Vector w = -wo;
        if (Vector.Dot(wo, n) < 0.0) {
          eta = 1.0 / eta;
        }

        Normal dndx = bsdf.dgShading.dndu * bsdf.dgShading.dudx +
                      bsdf.dgShading.dndv * bsdf.dgShading.dvdx;
        Normal dndy = bsdf.dgShading.dndu * bsdf.dgShading.dudy +
                      bsdf.dgShading.dndv * bsdf.dgShading.dvdy;

        Vector dwodx = -ray.rxDirection - wo;
        Vector dwody = -ray.ryDirection - wo;
        double dDNdx = Vector.Dot(dwodx, n) + Vector.Dot(wo, dndx);
        double dDNdy = Vector.Dot(dwody, n) + Vector.Dot(wo, dndy);

        double mu = eta * Vector.Dot(w, n) - Vector.Dot(wi, n);
        double dmudx = (eta - (eta * eta * Vector.Dot(w, n)) /
                       Vector.Dot(wi, n)) * dDNdx;
        double dmudy = (eta - (eta * eta * Vector.Dot(w, n)) /
                       Vector.Dot(wi, n)) * dDNdy;

        rd.rxDirection = wi + dwodx * eta -
                         new Vector.from(dndx * mu + n * dmudx);
        rd.ryDirection = wi + dwody * eta -
                         new Vector.from(dndy * mu + n * dmudy);
      }

      Stats.STARTED_SPECULAR_REFRACTION_RAY(rd);
      Spectrum Li = renderer.Li(scene, rd, sample, rng);
      L = f * Li * (Vector.AbsDot(wi, n) / pdf[0]);
      Stats.FINISHED_SPECULAR_REFRACTION_RAY(rd);
    }

    return L;
  }

  static Distribution1D ComputeLightSamplingCDF(Scene scene) {
    int nLights = scene.lights.length;

    Float32List lightPower = new Float32List(nLights);

    for (int i = 0; i < nLights; ++i) {
      lightPower[i] = scene.lights[i].power(scene).luminance();
    }

    return new Distribution1D(lightPower, nLights);
  }
}

