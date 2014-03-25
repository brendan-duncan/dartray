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
part of core;

/**
 * Defines the origin of a light source and the distribution of energy that
 * it emits.
 */
abstract class Light {
  Light(Transform l2w, [int ns = 1]) :
    nSamples = Math.max(1, ns),
    lightToWorld = new Transform.from(l2w),
    worldToLight = new Transform.from(Transform.Inverse(l2w)) {
      // Warn if light has transformation with scale
      if (worldToLight.hasScale()) {
        LogWarning("Scaling detected in world to light transformation!\n"
                   "The system has numerous assumptions, implicit and explicit,\n"
                   "that this transform will have no scale factors in it.\n"
                   "Proceed at your own risk; your image may have errors or\n"
                   "the system may crash as a result of this.");
      }
  }

  Spectrum sampleL(Point p, double pEpsilon,
        LightSample ls, double time, Vector wi, List<double> pdf,
        VisibilityTester vis);

  Spectrum power(Scene scene);

  bool isDeltaLight();

  Spectrum Le(RayDifferential r) {
    return new Spectrum(0.0);
  }

  double pdf(Point p, Vector wi);

  Spectrum sampleL2(Scene scene, LightSample ls,
                         double u1, double u2, double time,
                         Ray ray, Normal Ns, List<double> pdf);

  void shProject(Point p, double pEpsilon, int lmax,
                   Scene scene, bool computeLightVisibility, double time,
                   RNG rng, List<Spectrum> coeffs) {
    for (int i = 0; i < SHTerms(lmax); ++i) {
      coeffs[i] = new Spectrum(0.0);
    }

    int ns = RoundUpPow2(nSamples);
    int scramble1D = rng.randomUInt();
    List<int> scramble2D = [ rng.randomUInt(), rng.randomUInt() ];
    Float32List Ylm = new Float32List(SHTerms(lmax));

    List<double> u = [0.0, 0.0];
    List<double> pdf = [0.0];
    for (int i = 0; i < ns; ++i) {
      // Compute incident radiance sample from _light_, update SH _coeffs_
      Sample02(i, scramble2D, u);
      LightSample lightSample = new LightSample.set(u[0], u[1],
                                                  VanDerCorput(i, scramble1D));
      Vector wi = new Vector();
      VisibilityTester vis = new VisibilityTester();
      Spectrum Li = sampleL(p, pEpsilon, lightSample, time, wi, pdf, vis);
      if (!Li.isBlack() && pdf[0] > 0.0 &&
          (!computeLightVisibility || vis.unoccluded(scene))) {
        // Add light sample contribution to MC estimate of SH coefficients
        SHEvaluate(wi, lmax, Ylm);
        for (int j = 0; j < SHTerms(lmax); ++j) {
          coeffs[j] += Li * Ylm[j] / (pdf[0] * ns);
        }
      }
    }
  }

  final int nSamples;
  final Transform lightToWorld;
  final Transform worldToLight;
}
