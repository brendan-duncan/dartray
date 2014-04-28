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

class DistantLight extends Light {
  DistantLight(Transform light2world, Spectrum radiance, Vector dir)
      : super(light2world) {
    lightDir = Vector.Normalize(lightToWorld.transformVector(dir));
    L = radiance;
  }

  bool isDeltaLight() {
    return true;
  }

  Spectrum power(Scene scene) {
    Point worldCenter = new Point();
    double worldRadius = scene.worldBound.boundingSphere(worldCenter);
    return L * Math.PI * worldRadius * worldRadius;
  }

  Spectrum sampleLAtPoint(Point p, double pEpsilon, LightSample ls,
                          double time, Vector wi, List<double> pdf,
                          VisibilityTester vis) {
    wi.copy(lightDir);
    pdf[0] = 1.0;
    vis.setRay(p, pEpsilon, wi, time);
    return L;
  }

  Spectrum sampleL(Scene scene, LightSample ls, double u1, double u2,
                   double time, Ray ray, Normal Ns, List<double> pdf) {
    // Choose point on disk oriented toward infinite light direction
    Point worldCenter = new Point();
    double worldRadius = scene.worldBound.boundingSphere(worldCenter);

    Vector v1 = new Vector();
    Vector v2 = new Vector();
    Vector.CoordinateSystem(lightDir, v1, v2);

    List<double> d1 = [0.0];
    List<double> d2 = [0.0];
    ConcentricSampleDisk(ls.uPos[0], ls.uPos[1], d1, d2);
    Point Pdisk = worldCenter + (v1 * d1[0] + v2 * d2[0]) * worldRadius;

    // Set ray origin and direction for infinite light ray
    ray.set(Pdisk + lightDir * worldRadius, -lightDir, 0.0, INFINITY, time);
    Ns.copy(ray.direction);

    pdf[0] = 1.0 / (Math.PI * worldRadius * worldRadius);

    return L;
  }

  double pdf(Point p, Vector w) {
    return 0.0;
  }

  static DistantLight Create(Transform light2world, ParamSet paramSet) {
    Spectrum L = paramSet.findOneSpectrum('L', new Spectrum(1.0));
    Spectrum sc = paramSet.findOneSpectrum('scale', new Spectrum(1.0));
    Point from = paramSet.findOnePoint('from', new Point(0.0, 0.0, 0.0));
    Point to = paramSet.findOnePoint('to', new Point(0.0, 0.0, 1.0));
    Vector dir = from - to;
    return new DistantLight(light2world, L * sc, dir);
  }

  Vector lightDir;
  Spectrum L;
}
