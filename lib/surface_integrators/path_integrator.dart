/****************************************************************************
 *  Copyright (C) 2014 by Brendan Duncan.                                   *
 *                                                                          *
 *  This file is part of DartRay.                                           *
 *                                                                          *
 *  Licensed under the Apache License, Version 2.0 (the "License");         *
 *  you may not use this file except in compliance with the License.        *
 *  You may obtain a copy of the License at                                 *
 *                                                                          *
 *  http://www.apache.org/licenses/LICENSE-2.0                              *
 *                                                                          *
 *  Unless required by applicable law or agreed to in writing, software     *
 *  distributed under the License is distributed on an "AS IS" BASIS,       *
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 *  See the License for the specific language governing permissions and     *
 *  limitations under the License.                                          *
 *                                                                          *
 *  This project is based on PBRT v2 ; see http://www.pbrt.org              *
 *  pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys. *
 ****************************************************************************/
part of surface_integrators;

/**
 * Path tracing surface integrator.
 */
class PathIntegrator extends SurfaceIntegrator {
  PathIntegrator(int md) {
    maxDepth = md;
  }

  Spectrum Li(Scene scene, Renderer renderer,
              RayDifferential r, Intersection isect, Sample sample,
              RNG rng) {
    // Declare common path integration variables
    Spectrum pathThroughput = new Spectrum(1.0);
    Spectrum L = new Spectrum(0.0);

    RayDifferential ray = new RayDifferential.fromRay(r);
    bool specularBounce = false;
    Intersection localIsect = new Intersection();
    Intersection isectP = isect;

    Vector wi = new Vector();
    List<double> pdf = [0.0];
    List<int> flags = [0];

    for (int bounces = 0; ; ++bounces) {
      // Possibly add emitted light at path vertex
      if (bounces == 0 || specularBounce) {
        L += pathThroughput * isectP.Le(-ray.direction);
      }

      // Sample illumination from lights to find path contribution
      BSDF bsdf = isectP.getBSDF(ray);
      Point p = bsdf.dgShading.p;
      Normal n = bsdf.dgShading.nn;
      Vector wo = -ray.direction;

      if (bounces < SAMPLE_DEPTH) {
        L += pathThroughput *
            Integrator.UniformSampleOneLight(scene, renderer, p, n, wo,
                                   isectP.rayEpsilon, ray.time, bsdf,
                                   sample, rng,
                    lightNumOffset[bounces], lightSampleOffsets[bounces],
                    bsdfSampleOffsets[bounces]);
      } else {
        L += pathThroughput *
            Integrator.UniformSampleOneLight(scene, renderer, p, n, wo,
                                   isectP.rayEpsilon, ray.time, bsdf, sample,
                                   rng);
      }

      // Sample BSDF to get new path direction

      // Get _outgoingBSDFSample_ for sampling new path direction
      BSDFSample outgoingBSDFSample;
      if (bounces < SAMPLE_DEPTH) {
        outgoingBSDFSample = new BSDFSample.sample(sample,
                                                 pathSampleOffsets[bounces], 0);
      } else {
        outgoingBSDFSample = new BSDFSample.random(rng);
      }

      Spectrum f = bsdf.sample_f(wo, wi, outgoingBSDFSample, pdf, BSDF_ALL,
                                 flags);
      if (f.isBlack() || pdf[0] == 0.0) {
        break;
      }

      specularBounce = (flags[0] & BSDF_SPECULAR) != 0;
      pathThroughput *= f * Vector.AbsDot(wi, n) / pdf[0];

      ray = new RayDifferential.child(p, wi, ray, isectP.rayEpsilon);

      // Possibly terminate the path
      if (bounces > 3) {
        double continueProbability = Math.min(0.5, pathThroughput.luminance());
        if (rng.randomFloat() > continueProbability) {
          break;
        }
        pathThroughput /= continueProbability;
      }

      if (bounces == maxDepth) {
        break;
      }

      // Find next vertex of path
      if (!scene.intersect(ray, localIsect)) {
        if (specularBounce) {
          for (int i = 0; i < scene.lights.length; ++i) {
            L += pathThroughput * scene.lights[i].Le(ray);
          }
        }

        break;
      }

      pathThroughput *= renderer.transmittance(scene, ray, null, rng);

      isectP = localIsect;
    }

    return L;
  }

  void requestSamples(Sampler sampler, Sample sample, Scene scene) {
    for (int i = 0; i < SAMPLE_DEPTH; ++i) {
      lightSampleOffsets[i] = new LightSampleOffsets(1, sample);
      lightNumOffset[i] = sample.add1D(1);
      bsdfSampleOffsets[i] = new BSDFSampleOffsets(1, sample);
      pathSampleOffsets[i] = new BSDFSampleOffsets(1, sample);
    }
  }

  static PathIntegrator Create(ParamSet params) {
    int maxDepth = params.findOneInt('maxdepth', 5);
    return new PathIntegrator(maxDepth);
  }

  int maxDepth;
  static const int SAMPLE_DEPTH = 3;

  List<LightSampleOffsets> lightSampleOffsets =
      new List<LightSampleOffsets>(SAMPLE_DEPTH);

  List<int> lightNumOffset = new List<int>(SAMPLE_DEPTH);

  List<BSDFSampleOffsets> bsdfSampleOffsets =
      new List<BSDFSampleOffsets>(SAMPLE_DEPTH);

  List<BSDFSampleOffsets> pathSampleOffsets =
      new List<BSDFSampleOffsets>(SAMPLE_DEPTH);
}
