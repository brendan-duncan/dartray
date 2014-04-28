/****************************************************************************
 * Copyright (C) 2014 by Brendan Duncan.                                    *
 *                                                                          *
 * This file is part of DartRay.                                            *
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the 'License');          *
 * you may not use this file except in compliance with the License.         *
 * You may obtain a copy of the License at                                  *
 *                                                                          *
 * http://www.apache.org/licenses/LICENSE-2.0                               *
 *                                                                          *
 * Unless required by applicable law or agreed to in writing, software      *
 * distributed under the License is distributed on an 'AS IS' BASIS,        *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 * See the License for the specific language governing permissions and      *
 * limitations under the License.                                           *
 *                                                                          *
 * This project is based on PBRT v2 ; see http://www.pbrt.org               *
 * pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.  *
 ****************************************************************************/
part of core;

/**
 * Defines the origin of a light source and the distribution of energy that
 * it emits.
 */
abstract class Light {
  final Transform lightToWorld;
  final Transform worldToLight;
  /// The number of samples to evaluate on the surface of an area light.
  final int nSamples;

  Light(Transform l2w, [int ns = 1])
      : nSamples = Math.max(1, ns),
        lightToWorld = new Transform.from(l2w),
        worldToLight = new Transform.from(Transform.Inverse(l2w)) {
      // Warn if light has transformation with scale
      if (worldToLight.hasScale()) {
        LogWarning('Scaling detected in world to light transformation! '
                   'The system has numerous assumptions, implicit and '
                   'explicit, that this transform will have no scale factors '
                   'in it. Proceed at your own risk; your image may have '
                   'errors or the system may crash as a result of this.');
      }
  }

  /**
   * Returns the total emitted power of the light into the scene.
   * This quantity is useful for light transport algorithms that may want
   * to devote additional computational resources to lights in the scene
   * that make the largest distributions.
   */
  Spectrum power(Scene scene);

  /**
   * Indicates whether the light is described by a delta distribution.
   * Such lights include point lights, which emit illumination from a single
   * point, and directional lights, where all light arrives from the same
   * direction. The only way to detect illumination from light sources like
   * these is to call their sampleL methods. It is impossible to
   * randomly choose a direction from a point p that happens to find such
   * a light source.
   */
  bool isDeltaLight();

  /**
   * Returns the emitted radiance of the light along a ray that didn't hit
   * anything in the scene.
   */
  Spectrum Le(RayDifferential r) {
    return new Spectrum(0.0);
  }

  /**
   * The probability density function (PDF) describes the relative
   * probability of a random variable taking on a particular value.
   */
  double pdf(Point p, Vector wi);

  /**
   * Returns the radiance arriving at point [p] from the light.
   */
  Spectrum sampleLAtPoint(Point p, double pEpsilon, LightSample ls,
                          double time, Vector wi, List<double> pdf,
                          VisibilityTester vis);

  /**
   * Returns the radiance from the light in a direction ([u1],[u2])
   * from the sample point [ls] on the light.
   */
  Spectrum sampleL(Scene scene, LightSample ls,
                   double u1, double u2, double time,
                   Ray ray, Normal Ns, List<double> pdf);

  /**
   * Compute the spherical hamonic projection of the light.
   */
  void shProject(Point p, double pEpsilon, int lmax,
                 Scene scene, bool computeLightVisibility, double time,
                 RNG rng, List<Spectrum> coeffs) {
    for (int i = 0, len = SphericalHarmonics.Terms(lmax); i < len; ++i) {
      coeffs[i] = new Spectrum(0.0);
    }

    int ns = RoundUpPow2(nSamples);
    int scramble1D = rng.randomUint();
    List<int> scramble2D = [rng.randomUint(), rng.randomUint()];
    Float32List Ylm = new Float32List(SphericalHarmonics.Terms(lmax));

    List<double> u = [0.0, 0.0];
    List<double> pdf = [0.0];

    for (int i = 0; i < ns; ++i) {
      // Compute incident radiance sample from light, update SH coeffs
      Sample02(i, scramble2D, u);
      LightSample lightSample = new LightSample(u[0], u[1],
                                                VanDerCorput(i, scramble1D));
      Vector wi = new Vector();
      VisibilityTester vis = new VisibilityTester();
      Spectrum Li = sampleLAtPoint(p, pEpsilon, lightSample, time, wi, pdf,
                                   vis);
      if (!Li.isBlack() && pdf[0] > 0.0 &&
          (!computeLightVisibility || vis.unoccluded(scene))) {
        // Add light sample contribution to MC estimate of SH coefficients
        SphericalHarmonics.Evaluate(wi, lmax, Ylm);
        for (int j = 0; j < SphericalHarmonics.Terms(lmax); ++j) {
          coeffs[j] += Li * Ylm[j] / (pdf[0] * ns);
        }
      }
    }
  }
}
