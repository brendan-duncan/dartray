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

class Vector {
  double x, y, z;

  Vector([this.x = 0.0, this.y = 0.0, this.z = 0.0]);

  Vector.from(Vector other) :
    x = other.x,
    y = other.y,
    z = other.z;

  void copy(Vector other) {
    x = other.x;
    y = other.y;
    z = other.z;
  }

  Vector operator+(Vector v) =>
    new Vector(x + v.x, y + v.y, z + v.z);

  void add(Vector v) {
    x += v.x;
    y += v.y;
    z += v.z;
  }

  Vector operator-(Vector v) =>
      new Vector(x - v.x, y - v.y, z - v.z);

  void subtract(Vector v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
  }

  Vector operator*(double f) =>
      new Vector(x * f, y * f, z * f);

  void scale(double f) {
    x *= f;
    y *= f;
    z *= f;
  }

  Vector operator/(double f) =>
      new Vector(x / f, y / f, z / f);

  void invScale(double f) {
    x /= f;
    y /= f;
    z /= f;
  }

  Vector operator-() =>
      new Vector(-x, -y, -z);

  double operator[](int i) =>
      (i == 0) ? x : (i == 1) ? y : z;

  operator[]=(int i, double v) =>
    (i == 0) ? x = v : (i == 1) ? y = v : z = v;

  double lengthSquared() =>
      x * x + y * y + z * z;

  double length() =>
      Math.sqrt(lengthSquared());

  String toString() {
    return '$x $y $z';
  }

  static double CosTheta(Vector v) => v.z;

  static double AbsCosTheta(Vector v) => v.z.abs();

  static double SinTheta2(Vector v) =>
    Math.max(0.0, 1.0 - CosTheta(v) * CosTheta(v));

  static double SinTheta(Vector v) => Math.sqrt(SinTheta2(v));

  static double CosPhi(Vector v) {
    double sintheta = SinTheta(v);
    if (sintheta == 0.0) {
      return 1.0;
    }
    return (v.x / sintheta).clamp(-1.0, 1.0);
  }

  static double SinPhi(Vector v) {
    double sintheta = SinTheta(v);
    if (sintheta == 0.0) {
      return 0.0;
    }
    return (v.y / sintheta).clamp(-1.0, 1.0);
  }

  static double Distance(Vector a, Vector b) {
      return (b - a).length();
    }

  static double DistanceSquared(Vector a, Vector b) {
    return (b - a).lengthSquared();
  }

  static double Dot(Vector v1, Vector v2) =>
      v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;

  static double AbsDot(Vector v1, Vector v2) =>
      (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z).abs();

  static Vector Cross(Vector v1, Vector v2) {
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return new Vector((v1y * v2z) - (v1z * v2y),
                      (v1z * v2x) - (v1x * v2z),
                      (v1x * v2y) - (v1y * v2x));
  }

  static Vector Normalize(Vector v) => v / v.length();

  static Vector SphericalDirection(double sintheta, double costheta,
                                   double phi) {
      return new Vector(sintheta * Math.cos(phi),
                       sintheta * Math.sin(phi),
                       costheta);
  }

  static Vector SphericalDirectionVec(double sintheta, double costheta,
                                      double phi,
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
      double invLen = 1.0 / Math.sqrt(v1.x * v1.x + v1.z * v1.z);
      v2.x = -v1.z * invLen;
      v2.y = 0.0;
      v2.z = v1.x * invLen;
    } else {
      double invLen = 1.0 / Math.sqrt(v1.y * v1.y + v1.z * v1.z);
      v2.x = 0.0;
      v2.y = v1.z * invLen;
      v2.z = -v1.y * invLen;
    }

    v3.copy(Vector.Cross(v1, v2));
  }

  static Vector FaceForward(Vector n, Vector n2) {
    return (Vector.Dot(n, n2) < 0.0) ? new Vector(-n.x, -n.y, -n.z) : n;
  }
}
