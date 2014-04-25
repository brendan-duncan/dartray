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

class DiffusePRTIntegrator extends SurfaceIntegrator {
  DiffusePRTIntegrator(int lm, int ns) :
    lmax = lm,
    nSamples = RoundUpPow2(ns),
    c_in = new List<Spectrum>(SphericalHarmonics.Terms(lm)) {
    for (int i = 0, len = c_in.length; i < len; ++i) {
      c_in[i] = new Spectrum(0.0);
    }
  }

  void preprocess(Scene scene, Camera camera, Renderer renderer) {
    BBox bbox = scene.worldBound;
    Point p = bbox.pMin * 0.5 + bbox.pMax * 0.5;
    RNG rng = new RNG();
    SphericalHarmonics.ProjectIncidentDirectRadiance(p, 0.0, camera.shutterOpen,
                                                     scene, false, lmax, rng,
                                                     c_in);
  }

  void requestSamples(Sampler sampler, Sample sample, Scene scene) {
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
    // Compute reflected radiance using diffuse PRT

    // Project diffuse transfer function at point to SH
    List<Spectrum> c_transfer = new List<Spectrum>(SphericalHarmonics.Terms(lmax));
    for (int i = 0, len = c_transfer.length; i < len; ++i) {
      c_transfer[i] = new Spectrum(0.0);
    }

    SphericalHarmonics.ComputeDiffuseTransfer(p, Normal.FaceForward(n, wo),
                                              isect.rayEpsilon,
                                              scene, rng, nSamples, lmax,
                                              c_transfer);

    // Compute integral of product of incident radiance and transfer function
    Spectrum Kd = bsdf.rho2(wo, rng, BSDF_ALL_REFLECTION) * INV_PI;
    Spectrum Lo = new Spectrum(0.0);
    for (int i = 0, len = SphericalHarmonics.Terms(lmax); i < len; ++i) {
      Lo += c_in[i] * c_transfer[i];
    }

    return L + Kd * Lo.clamp();
  }

  static DiffusePRTIntegrator Create(ParamSet params) {
    int lmax = params.findOneInt("lmax", 4);
    int ns = params.findOneInt("nsamples", 4096);
    return new DiffusePRTIntegrator(lmax, ns);
  }

  final int lmax;
  final int nSamples;
  final List<Spectrum> c_in;
}
