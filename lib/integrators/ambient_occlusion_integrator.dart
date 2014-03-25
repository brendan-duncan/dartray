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
 *   This project is based on PBRT v2 ; see http://www.pbrt.org             *
 *   pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.*
 ****************************************************************************/
part of integrators;

class AmbientOcclusionIntegrator extends SurfaceIntegrator {
  AmbientOcclusionIntegrator(int ns, this.maxDist) {
    nSamples = RoundUpPow2(ns);
  }

  static AmbientOcclusionIntegrator Create(ParamSet params) {
    int nSamples = params.findOneInt("nsamples", 2048);
    double maxDist = params.findOneFloat("maxdist", INFINITY);
    return new AmbientOcclusionIntegrator(nSamples, maxDist);
  }

  Spectrum Li(Scene scene, Renderer renderer,
                 RayDifferential ray, Intersection isect,
                 Sample sample, RNG rng) {
    BSDF bsdf = isect.getBSDF(ray);
    Point p = bsdf.dgShading.p;
    Normal n = Normal.FaceForward(isect.dg.nn, -ray.direction);

    List<int> scramble = [ rng.randomUInt(), rng.randomUInt() ];
    List<double> u = [0.0, 0.0];
    int nClear = 0;
    for (int i = 0; i < nSamples; ++i) {
      Sample02(i, scramble, u);
      Vector w = UniformSampleSphere(u[0], u[1]);
      if (Vector.Dot(w, n) < 0.0) {
        w = -w;
      }

      Ray r = new Ray(p, w, 0.01, maxDist);

      if (!scene.intersectP(r)) {
        ++nClear;
      }
    }

    return new Spectrum(nClear / nSamples);
  }

  int nSamples;
  double maxDist;
}
