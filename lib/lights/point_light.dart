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
part of lights;

class PointLight extends Light {
  PointLight(Transform light2world, this.intensity)
      : super(light2world) {
    lightPos = light2world.transformPoint(new Point());
  }

  Spectrum sampleLAtPoint(Point p, double pEpsilon, LightSample ls,
        double time, Vector wi, List<double> pdf, VisibilityTester vis) {
    wi.copy(Vector.Normalize(lightPos - p));
    pdf[0] = 1.0;
    vis.setSegment(p, pEpsilon, lightPos, 0.0, time);
    return intensity / Vector.DistanceSquared(lightPos, p);
  }

  Spectrum power(Scene scene) {
    return intensity * (4.0 * Math.PI);
  }

  bool isDeltaLight() {
    return true;
  }

  Spectrum sampleL(Scene scene, LightSample ls, double u1,
                   double u2, double time, Ray ray, Normal Ns,
                   List<double> pdf) {
    ray.origin = new Point.from(lightPos);
    ray.direction = UniformSampleSphere(ls.uPos[0], ls.uPos[1]);
    ray.minDistance = 0.0;
    ray.maxDistance = INFINITY;
    ray.time = time;

    Ns.copy(ray.direction);
    pdf[0] = UniformSpherePdf();
    return new Spectrum.from(intensity);
  }

  double pdf(Point p, Vector w) {
    return 0.0;
  }

  void shProject(Point p, double pEpsilon, int lmax, Scene scene,
        bool computeLightVisibility, double time, RNG rng,
        List<Spectrum> coeffs) {
    for (int i = 0; i < SphericalHarmonics.Terms(lmax); ++i) {
      coeffs[i] = new Spectrum(0.0);
    }

    if (computeLightVisibility &&
        scene.intersectP(new Ray(p, Vector.Normalize(lightPos - p),
                                 pEpsilon, Vector.Distance(lightPos, p),
                                 time))) {
      return;
    }

    // Project point light source to SH
    List<double> Ylm = new List<double>(SphericalHarmonics.Terms(lmax));
    Vector wi = Vector.Normalize(lightPos - p);
    SphericalHarmonics.Evaluate(wi, lmax, Ylm);
    Spectrum Li = intensity / Vector.DistanceSquared(lightPos, p);
    for (int i = 0; i < SphericalHarmonics.Terms(lmax); ++i) {
      coeffs[i] = Li * Ylm[i];
    }
  }

  static PointLight Create(Transform light2world, ParamSet paramSet) {
    Spectrum I = paramSet.findOneSpectrum('I', new Spectrum(1.0));
    Spectrum sc = paramSet.findOneSpectrum('scale', new Spectrum(1.0));
    Point P = paramSet.findOnePoint('from', new Point());
    Transform l2w = Transform.Translate(P) * light2world;
    return new PointLight(l2w, I * sc);
  }

  Point lightPos;
  Spectrum intensity;
}
