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
part of core;

/**
 * A 3-dimensional vector.
 */
class Vector {
  final Float32List data;

  Vector([num x = 0.0, num y = 0.0, num z = 0.0])
      : data = new Float32List(3) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
  }

  Vector.from(Vector other)
      : data = new Float32List.fromList(other.data);

  double get x => data[0];

  set x(num v) => data[0] = v;

  double get y => data[1];

  set y(num v) => data[1] = v;

  double get z => data[2];

  set z(num v) => data[2] = v;

  void copy(Vector other) {
    data[0] = other.data[0];
    data[1] = other.data[1];
    data[2] = other.data[2];
  }

  Vector operator +(Vector v) =>
      new Vector(data[0] + v.data[0],
                 data[1] + v.data[1],
                 data[2] + v.data[2]);

  Vector operator -(Vector v) =>
      new Vector(data[0] - v.data[0],
                 data[1] - v.data[1],
                 data[2] - v.data[2]);

  Vector operator *(num f) =>
      new Vector(data[0] * f, data[1] * f, data[2] * f);

  Vector operator /(num f) =>
      new Vector(data[0] / f, data[1] / f, data[2] / f);

  Vector operator -() =>
      new Vector(-data[0], -data[1], -data[2]);

  double operator [](int i) => data[i];

  operator []=(int i, num v) => data[i] = v;

  double lengthSquared() =>
      data[0] * data[0] + data[1] * data[1] + data[2] * data[2];

  double length() => Math.sqrt(lengthSquared());

  Vector scale(num s) {
    data[0] *= s;
    data[1] *= s;
    data[2] *= s;
    return this;
  }

  Vector invScale(num s) {
    data[0] /= s;
    data[1] /= s;
    data[2] /= s;
    return this;
  }

  Vector add(Vector s) {
    data[0] += s.data[0];
    data[1] += s.data[1];
    data[2] += s.data[2];
    return this;
  }

  Vector subtract(Vector s) {
    data[0] -= s.data[0];
    data[1] -= s.data[1];
    data[2] -= s.data[2];
    return this;
  }

  String toString() => '${data[0]} ${data[1]} ${data[2]}';

  bool hasNaNs() => data[0].isNaN || data[1].isNaN || data[2].isNaN;

  static double CosTheta(Vector v) => v.data[2];

  static double AbsCosTheta(Vector v) => v.data[2].abs();

  static double SinTheta2(Vector v) =>
    Math.max(0.0, 1.0 - CosTheta(v) * CosTheta(v));

  static double SinTheta(Vector v) => Math.sqrt(SinTheta2(v));

  static double CosPhi(Vector v) {
    double sintheta = SinTheta(v);
    if (sintheta == 0.0) {
      return 1.0;
    }
    return (v.data[0] / sintheta).clamp(-1.0, 1.0);
  }

  static double SinPhi(Vector v) {
    double sintheta = SinTheta(v);
    if (sintheta == 0.0) {
      return 0.0;
    }
    return (v.data[1] / sintheta).clamp(-1.0, 1.0);
  }

  static double Distance(Vector a, Vector b) {
    return (b - a).length();
  }

  static double DistanceSquared(Vector a, Vector b) {
    return (b - a).lengthSquared();
  }

  static double Dot(Vector v1, Vector v2) =>
      v1.data[0] * v2.data[0] + v1.data[1] * v2.data[1] +
      v1.data[2] * v2.data[2];

  static double AbsDot(Vector v1, Vector v2) =>
      (v1.data[0] * v2.data[0] + v1.data[1] * v2.data[1] +
       v1.data[2] * v2.data[2]).abs();

  static Vector Cross(Vector v1, Vector v2) {
    double v1x = v1.data[0];
    double v1y = v1.data[1];
    double v1z = v1.data[2];
    double v2x = v2.data[0];
    double v2y = v2.data[1];
    double v2z = v2.data[2];
    return new Vector((v1y * v2z) - (v1z * v2y),
                      (v1z * v2x) - (v1x * v2z),
                      (v1x * v2y) - (v1y * v2x));
  }

  static Vector Normalize(Vector v) => v / v.length();

  static Vector SphericalDirection(num sintheta, num costheta, num phi) {
      return new Vector(sintheta * Math.cos(phi),
                       sintheta * Math.sin(phi),
                       costheta);
  }

  static Vector SphericalDirectionVec(num sintheta, num costheta, num phi,
                                      Vector x, Vector y, Vector z) {
      return x * (sintheta * Math.cos(phi)) +
             y * (sintheta * Math.sin(phi)) +
             z * costheta;
  }

  static double SphericalTheta(Vector v) {
    return Math.acos(v.z.clamp(-1.0, 1.0));
  }

  static double SphericalPhi(Vector v) {
    double p = Math.atan2(v.y, v.x);
    return (p < 0.0) ? p + 2.0 * Math.PI : p;
  }

  static bool SameHemisphere(Vector w, Vector wp) {
    return w.z * wp.z > 0.0;
  }

  static void CoordinateSystem(Vector v1, Vector v2, Vector v3) {
    if (v1.x.abs() > v1.y.abs()) {
      double invLen = 1.0 / Math.sqrt(v1.data[0] * v1.data[0] +
                                      v1.data[2] * v1.data[2]);
      v2.data[0] = -v1.data[2] * invLen;
      v2.data[1] = 0.0;
      v2.data[2] = v1.data[0] * invLen;
    } else {
      double invLen = 1.0 / Math.sqrt(v1.data[1] * v1.data[1] +
                                      v1.data[2] * v1.data[2]);
      v2.data[0] = 0.0;
      v2.data[1] = v1.data[2] * invLen;
      v2.data[2] = -v1.data[1] * invLen;
    }

    v3.copy(Vector.Cross(v1, v2));
  }

  static Vector FaceForward(Vector n, Vector n2) {
    return (Vector.Dot(n, n2) < 0.0) ?
           new Vector(-n.data[0], -n.data[1], -n.data[2]) :
           n;
  }
}
