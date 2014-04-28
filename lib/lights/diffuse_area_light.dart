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

class DiffuseAreaLight extends AreaLight {
  DiffuseAreaLight(Transform light2world,
                   Spectrum Le, int ns, Shape shape)
      : Lemit = new Spectrum.from(Le),
        shapeSet = new ShapeSet(shape),
        super(light2world, ns) {
    area = shapeSet.area;
  }

  Spectrum L(Point p, Normal n, Vector w) {
    return Vector.Dot(n, w) > 0.0 ? Lemit : new Spectrum(0.0);
  }

  Spectrum power(Scene scene) {
    return Lemit * area * Math.PI;
  }

  bool isDeltaLight() {
    return false;
  }

  double pdf(Point p, Vector w) {
    return shapeSet.pdf(p, w);
  }

  Spectrum sampleLAtPoint(Point p, double pEpsilon, LightSample ls, double time,
                          Vector wo, List<double> pdf,
                          VisibilityTester visibility) {
    Normal ns = new Normal();
    Point ps = shapeSet.sample(ls, ns, p);
    wo.copy(Vector.Normalize(ps - p));
    pdf[0] = shapeSet.pdf(p, wo);
    visibility.setSegment(p, pEpsilon, ps, 1.0e-3, time);
    Spectrum Ls = L(ps, ns, -wo);
    return Ls;
  }

  Spectrum sampleL(Scene scene, LightSample ls, double u1, double u2,
                   double time, Ray ray, Normal Ns, List<double> pdf) {
    Point org = shapeSet.sample(ls, Ns);
    Vector dir = UniformSampleSphere(u1, u2);
    if (Vector.Dot(dir, Ns) < 0.0) {
      dir *= -1.0;
    }

    ray.origin = org;
    ray.direction = dir;
    ray.minDistance = 1.0e-3;
    ray.maxDistance = INFINITY;
    ray.time = time;
    pdf[0] = shapeSet.pdf(org) * INV_TWOPI;

    Spectrum Ls = L(org, Ns, dir);
    return Ls;
  }

  static DiffuseAreaLight Create(Transform light2world, ParamSet paramSet,
                                 Shape shape) {
    Spectrum L = paramSet.findOneSpectrum('L', new Spectrum(1.0));
    Spectrum sc = paramSet.findOneSpectrum('scale', new Spectrum(1.0));
    int nSamples = paramSet.findOneInt('nsamples', 1);
    return new DiffuseAreaLight(light2world, L * sc, nSamples, shape);
  }

  Spectrum Lemit;
  ShapeSet shapeSet;
  double area;
}
