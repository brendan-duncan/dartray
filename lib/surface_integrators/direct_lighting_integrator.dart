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

class DirectLightingIntegrator extends SurfaceIntegrator {
  static const int SAMPLE_ALL_UNIFORM = 0;
  static const int SAMPLE_ONE_UNIFORM = 1;

  DirectLightingIntegrator([this.strategy = SAMPLE_ALL_UNIFORM,
                            this.maxDepth = 5]);

  Spectrum Li(Scene scene, Renderer renderer,
      RayDifferential ray, Intersection isect,
      Sample sample, RNG rng) {
    Spectrum L = new Spectrum(0.0);
    // Evaluate BSDF at hit point
    BSDF bsdf = isect.getBSDF(ray);
    Vector wo = -ray.direction;
    Point p = bsdf.dgShading.p;
    Normal n = bsdf.dgShading.nn;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Compute direct lighting for _DirectLightingIntegrator_ integrator
    if (scene.lights.length > 0) {
      // Apply direct lighting strategy
      switch (strategy) {
        case SAMPLE_ALL_UNIFORM:
          L += Integrator.UniformSampleAllLights(scene, renderer, p, n, wo,
              isect.rayEpsilon, ray.time, bsdf, sample, rng,
              lightSampleOffsets, bsdfSampleOffsets);
          break;
        case SAMPLE_ONE_UNIFORM:
          L += Integrator.UniformSampleOneLight(scene, renderer, p, n, wo,
              isect.rayEpsilon, ray.time, bsdf, sample, rng,
              lightNumOffset, lightSampleOffsets[0], bsdfSampleOffsets[0]);
          break;
      }
    }

    if (ray.depth + 1 < maxDepth) {
      // Trace rays for specular reflection and refraction
      L += Integrator.SpecularReflect(ray, bsdf, rng, isect, renderer, scene,
                                      sample);
      L += Integrator.SpecularTransmit(ray, bsdf, rng, isect, renderer, scene,
                                       sample);
    }

    return L;
  }

  void requestSamples(Sampler sampler, Sample sample, Scene scene) {
    if (strategy == SAMPLE_ALL_UNIFORM) {
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
      lightNumOffset = -1;
    } else {
      // Allocate and request samples for sampling one light
      lightSampleOffsets = new List<LightSampleOffsets>(1);
      lightSampleOffsets[0] = new LightSampleOffsets(1, sample);
      lightNumOffset = sample.add1D(1);
      bsdfSampleOffsets = new List<BSDFSampleOffsets>(1);
      bsdfSampleOffsets[0] = new BSDFSampleOffsets(1, sample);
    }
  }

  static DirectLightingIntegrator Create(ParamSet params) {
    int maxDepth = params.findOneInt('maxdepth', 5);
    int strategy;
    String st = params.findOneString('strategy', 'all');
    if (st == 'one') {
      strategy = SAMPLE_ONE_UNIFORM;
    } else if (st == 'all') {
      strategy = SAMPLE_ALL_UNIFORM;
    } else {
      LogWarning('Strategy \'$st\' for direct lighting unknown. Using \'all\'.');
      strategy = SAMPLE_ALL_UNIFORM;
    }
    return new DirectLightingIntegrator(strategy, maxDepth);
  }

  int strategy;
  int maxDepth;
  List<LightSampleOffsets> lightSampleOffsets;
  List<BSDFSampleOffsets> bsdfSampleOffsets;
  int lightNumOffset;
}
