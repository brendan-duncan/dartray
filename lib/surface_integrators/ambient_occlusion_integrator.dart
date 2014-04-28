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

class AmbientOcclusionIntegrator extends SurfaceIntegrator {
  AmbientOcclusionIntegrator(int ns, this.minDist, this.maxDist) {
    nSamples = RoundUpPow2(ns);
  }

  Spectrum Li(Scene scene, Renderer renderer,
                 RayDifferential ray, Intersection isect,
                 Sample sample, RNG rng) {
    BSDF bsdf = isect.getBSDF(ray);
    Point p = bsdf.dgShading.p;
    Normal n = Normal.FaceForward(isect.dg.nn, -ray.direction);

    List<int> scramble = [rng.randomUint(), rng.randomUint()];
    List<double> u = [0.0, 0.0];
    int nClear = 0;
    for (int i = 0; i < nSamples; ++i) {
      Sample02(i, scramble, u);
      Vector w = UniformSampleSphere(u[0], u[1]);
      if (Vector.Dot(w, n) < 0.0) {
        w = -w;
      }

      Ray r = new Ray(p, w, minDist, maxDist);

      if (!scene.intersectP(r)) {
        ++nClear;
      }
    }

    return new Spectrum(nClear / nSamples);
  }

  static AmbientOcclusionIntegrator Create(ParamSet params) {
    int nSamples = params.findOneInt('nsamples', 2048);
    double maxDist = params.findOneFloat('maxdist', INFINITY);
    double minDist = params.findOneFloat('mindist', 1.0e-4);
    return new AmbientOcclusionIntegrator(nSamples, minDist, maxDist);
  }

  int nSamples;
  double minDist;
  double maxDist;
}
