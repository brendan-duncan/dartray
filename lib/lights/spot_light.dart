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

class SpotLight extends Light {
  SpotLight(Transform light2world, Spectrum I, double width, double fall)
      : super(light2world) {
    lightPos = lightToWorld.transformPoint(new Point(0.0, 0.0, 0.0));
    intensity = I;
    cosTotalWidth = Math.cos(Radians(width));
    cosFalloffStart = Math.cos(Radians(fall));
  }

  bool isDeltaLight() {
    return true;
  }

  double falloff(Vector w) {
    Vector wl = Vector.Normalize(worldToLight.transformVector(w));
    double costheta = wl.z;
    if (costheta < cosTotalWidth) {
      return 0.0;
    }

    if (costheta > cosFalloffStart) {
      return 1.0;
    }

    // Compute falloff inside spotlight cone
    double delta = (costheta - cosTotalWidth) /
                   (cosFalloffStart - cosTotalWidth);

    return delta * delta * delta * delta;
  }

  Spectrum power(Scene s) {
    return intensity * 2.0 * Math.PI *
               (1.0 - 0.5 * (cosFalloffStart + cosTotalWidth));
  }

  Spectrum sampleLAtPoint(Point p, double pEpsilon, LightSample ls,
                          double time, Vector wi, List<double> pdf,
                          VisibilityTester vis) {
    wi.copy(Vector.Normalize(lightPos - p));
    pdf[0] = 1.0;
    vis.setSegment(p, pEpsilon, lightPos, 0.0, time);

    return intensity * falloff(-wi) / Vector.DistanceSquared(lightPos, p);
  }

  Spectrum sampleL(Scene scene, LightSample ls, double u1, double u2,
                   double time, Ray ray, Normal Ns, List<double> pdf) {
    Vector v = UniformSampleCone(ls.uPos[0], ls.uPos[1], cosTotalWidth);
    ray.set(lightPos, lightToWorld.transformVector(v), 0.0, INFINITY, time);
    Ns.copy(ray.direction);
    pdf[0] = UniformConePdf(cosTotalWidth);
    return intensity * falloff(ray.direction);
  }

  double pdf(Point p, Vector w) {
    return 0.0;
  }

  static SpotLight Create(Transform l2w, ParamSet paramSet) {
    Spectrum I = paramSet.findOneSpectrum('I', new Spectrum(1.0));
    Spectrum sc = paramSet.findOneSpectrum('scale', new Spectrum(1.0));
    double coneangle = paramSet.findOneFloat('coneangle', 30.0);
    double conedelta = paramSet.findOneFloat('conedeltaangle', 5.0);
    // Compute spotlight world to light transformation
    Point from = paramSet.findOnePoint('from', new Point(0.0, 0.0, 0.0));
    Point to = paramSet.findOnePoint('to', new Point(0.0, 0.0, 1.0));
    Vector dir = Vector.Normalize(to - from);
    Vector du = new Vector();
    Vector dv = new Vector();
    Vector.CoordinateSystem(dir, du, dv);

    Transform dirToZ =
        new Transform(new Matrix4x4.values(du.x,  du.y,  du.z,  0.0,
                                           dv.x,  dv.y,  dv.z,  0.0,
                                           dir.x, dir.y, dir.z, 0.0,
                                           0.0,   0.0,   0.0,   1.0));

    Transform light2world = l2w *
        Transform.Translate(new Vector(from.x, from.y, from.z)) *
            Transform.Inverse(dirToZ);

    return new SpotLight(light2world, I * sc, coneangle,
                         coneangle - conedelta);
  }

  Point lightPos;
  Spectrum intensity;
  double cosTotalWidth;
  double cosFalloffStart;
}
