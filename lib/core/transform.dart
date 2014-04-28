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
 * Defines the 3d transformation of an object, stores the transformation matrix
 * and it's inverse.
 */
class Transform {
  Matrix4x4 m;
  Matrix4x4 mInv;

  Transform([Matrix4x4 m, Matrix4x4 inv])
      : this.m = (m == null) ? new Matrix4x4() : new Matrix4x4.from(m),
        this.mInv = (inv == null) ?
                    (m == null ? new Matrix4x4() : Matrix4x4.Inverse(m)) :
                    new Matrix4x4.from(inv);

  Transform.from(Transform t)
      : m = new Matrix4x4.from(t.m),
        mInv = new Matrix4x4.from(t.mInv);

  Transform copy(Transform t) {
    m = new Matrix4x4.from(t.m);
    mInv = new Matrix4x4.from(t.mInv);
    return this;
  }

  bool isIdentity() {
    return (m.data[0] == 1.0 && m.data[1] == 0.0 &&
            m.data[2] == 0.0 && m.data[3] == 0.0 &&
            m.data[4] == 0.0 && m.data[5] == 1.0 &&
            m.data[6] == 0.0 && m.data[7] == 0.0 &&
            m.data[8] == 0.0 && m.data[9] == 0.0 &&
            m.data[10] == 1.0 && m.data[11] == 0.0 &&
            m.data[12] == 0.0 && m.data[13] == 0.0 &&
            m.data[14] == 0.0 && m.data[15] == 1.0);
  }

  static Transform Inverse(Transform t) {
    return new Transform(t.mInv, t.m);
  }

  static Transform Transpose(Transform t) {
    return new Transform(Matrix4x4.Transpose(t.m),
                         Matrix4x4.Transpose(t.mInv));
  }

  bool operator ==(Transform t) {
    return t.m == m && t.mInv == mInv;
  }

  bool operator <(Transform t2) {
    for (int i = 0; i < 16; ++i) {
      if (m.data[i] < t2.m.data[i]) {
        return true;
      }
      if (m.data[i] > t2.m.data[i]) {
        return false;
      }
    }
    return false;
  }

  Transform operator *(Transform t2) {
    return new Transform(Matrix4x4.Mul(m, t2.m),
                         Matrix4x4.Mul(t2.mInv, mInv));
  }

  bool hasScale() {
    double la2 = transformVector(new Vector(1.0, 0.0, 0.0)).lengthSquared();
    double lb2 = transformVector(new Vector(0.0, 1.0, 0.0)).lengthSquared();
    double lc2 = transformVector(new Vector(0.0, 0.0, 1.0)).lengthSquared();
    return (la2 < 0.999 || la2 > 1.001) ||
           (lb2 < 0.999 || lb2 > 1.001) ||
           (lc2 < 0.999 || lc2 > 1.001);
  }

  bool swapsHandedness() {
    double det = ((m.data[0] *
                   (m.data[5] * m.data[10] -
                    m.data[6] * m.data[5])) -
                   (m.data[1] *
                    (m.data[4] * m.data[10] -
                     m.data[6] * m.data[4])) +
                    (m.data[2] *
                     (m.data[4] * m.data[9] -
                      m.data[5] * m.data[8])));
    return det < 0.0;
  }

  Point transformPoint(Point p, [Point out]) {
    if (out == null) {
      out = new Point();
    }

    double x = p.x;
    double y = p.y;
    double z = p.z;

    out.x = m.data[0] * x + m.data[1] * y + m.data[2] * z + m.data[3];
    out.y = m.data[4] * x + m.data[5] * y + m.data[6] * z + m.data[7];
    out.z = m.data[8] * x + m.data[9] * y + m.data[10] * z + m.data[11];
    double w = m.data[12] * x + m.data[13] * y + m.data[14] * z + m.data[15];

    if (w != 1.0) {
      out.invScale(w);
    }

    return out;
  }

  Vector transformVector(Vector p, [Vector out]) {
    if (out == null) {
      out = new Vector();
    }

    double x = p.x;
    double y = p.y;
    double z = p.z;

    out.x = m.data[0] * x + m.data[1] * y + m.data[2] * z;
    out.y = m.data[4] * x + m.data[5] * y + m.data[6] * z;
    out.z = m.data[8] * x + m.data[9] * y + m.data[10] * z;

    return out;
  }

  Normal transformNormal(Normal p, [Normal out]) {
    if (out == null) {
      out = new Normal();
    }

    double x = p.x;
    double y = p.y;
    double z = p.z;

    out.x = mInv.data[0] * x + mInv.data[4] * y + mInv.data[8] * z;
    out.y = mInv.data[1] * x + mInv.data[5] * y + mInv.data[9] * z;
    out.z = mInv.data[2] * x + mInv.data[6] * y + mInv.data[10] * z;

    return out;
  }

  BBox transformBBox(BBox b, [BBox out]) {
    if (out == null) {
      out = new BBox();
    }

    out.setPoint(transformPoint(b.pMin));
    out.unionPoint(transformPoint(new Point(b.pMax.x, b.pMin.y, b.pMin.z)));
    out.unionPoint(transformPoint(new Point(b.pMin.x, b.pMax.y, b.pMin.z)));
    out.unionPoint(transformPoint(new Point(b.pMin.x, b.pMin.y, b.pMax.z)));
    out.unionPoint(transformPoint(new Point(b.pMin.x, b.pMax.y, b.pMax.z)));
    out.unionPoint(transformPoint(new Point(b.pMax.x, b.pMax.y, b.pMin.z)));
    out.unionPoint(transformPoint(new Point(b.pMax.x, b.pMin.y, b.pMax.z)));
    out.unionPoint(transformPoint(b.pMax));

    return out;
  }

  Ray transformRay(Ray r, [Ray tr]) {
    if (tr == null) {
      tr = new Ray();
    }

    transformPoint(r.origin, tr.origin);
    transformVector(r.direction, tr.direction);

    if (tr != r) {
      tr.minDistance = r.minDistance;
      tr.maxDistance = r.maxDistance;
      tr.time = r.time;
      tr.depth = r.depth;
    }

    return tr;
  }

  RayDifferential transformRayDifferential(RayDifferential r,
                                           [RayDifferential rt]) {
    if (rt == null) {
      rt = new RayDifferential();
    }

    transformRay(r, rt);
    rt.hasDifferentials = r.hasDifferentials;
    transformPoint(r.rxOrigin, rt.rxOrigin);
    transformPoint(r.ryOrigin, rt.ryOrigin);
    transformVector(r.rxDirection, rt.rxDirection);
    transformVector(r.ryDirection, rt.ryDirection);

    return rt;
  }

  static Transform Translate(Vector delta) {
    Matrix4x4 m = new Matrix4x4.values(
                  1.0, 0.0, 0.0, delta.x,
                  0.0, 1.0, 0.0, delta.y,
                  0.0, 0.0, 1.0, delta.z,
                  0.0, 0.0, 0.0, 1.0);
    Matrix4x4 minv = new Matrix4x4.values(
                   1.0, 0.0, 0.0, -delta.x,
                   0.0, 1.0, 0.0, -delta.y,
                   0.0, 0.0, 1.0, -delta.z,
                   0.0, 0.0, 0.0, 1.0);

    return new Transform(m, minv);
  }

  static Transform Scale(double x, double y, double z) {
    Matrix4x4 m = new Matrix4x4.values(
                x, 0.0, 0.0, 0.0,
                0.0, y, 0.0, 0.0,
                0.0, 0.0, z, 0.0,
                0.0, 0.0, 0.0, 1.0);
    Matrix4x4 minv = new Matrix4x4.values(
                1.0 / x,     0.0,     0.0,     0.0,
                0.0,     1.0 / y,     0.0,     0.0,
                0.0,         0.0,     1.0 / z, 0.0,
                0.0,         0.0,     0.0,     1.0);
    return new Transform(m, minv);
  }

  static Transform RotateX(double angle) {
    double sin_t = Math.sin(Radians(angle));
    double cos_t = Math.cos(Radians(angle));
    Matrix4x4 m = new Matrix4x4.values(
        1.0, 0.0,   0.0,    0.0,
        0.0, cos_t, -sin_t, 0.0,
        0.0, sin_t, cos_t,  0.0,
        0.0, 0.0,   0.0,    1.0);
    return new Transform(m, Matrix4x4.Transpose(m));
  }

  static Transform RotateY(double angle) {
    double sin_t = Math.sin(Radians(angle));
    double cos_t = Math.cos(Radians(angle));
    Matrix4x4 m = new Matrix4x4.values(
        cos_t,  0.0,  sin_t, 0.0,
        0.0,    1.0,  0.0,   0.0,
        -sin_t, 0.0,  cos_t, 0.0,
        0.0,    0.0,  0.0,   1.0);
    return new Transform(m, Matrix4x4.Transpose(m));
  }

  static Transform RotateZ(double angle) {
    double sin_t = Math.sin(Radians(angle));
    double cos_t = Math.cos(Radians(angle));
    Matrix4x4 m = new Matrix4x4.values(
                cos_t, -sin_t, 0.0, 0.0,
                sin_t,  cos_t, 0.0, 0.0,
                0.0,      0.0, 1.0, 0.0,
                0.0,      0.0, 0.0, 1.0);
    return new Transform(m, Matrix4x4.Transpose(m));
  }

  static Transform Rotate(double angle, Vector axis) {
    Vector a = Vector.Normalize(axis);
    double s = Math.sin(Radians(angle));
    double c = Math.cos(Radians(angle));
    Matrix4x4 m = new Matrix4x4();

    m[0] = a.x * a.x + (1.0 - a.x * a.x) * c;
    m[1] = a.x * a.y * (1.0 - c) - a.z * s;
    m[2] = a.x * a.z * (1.0 - c) + a.y * s;
    m[3] = 0.0;

    m[4] = a.x * a.y * (1.0 - c) + a.z * s;
    m[5] = a.y * a.y + (1.0 - a.y * a.y) * c;
    m[6] = a.y * a.z * (1.0 - c) - a.x * s;
    m[7] = 0.0;

    m[8] = a.x * a.z * (1.0 - c) - a.y * s;
    m[9] = a.y * a.z * (1.0 - c) + a.x * s;
    m[10] = a.z * a.z + (1.0 - a.z * a.z) * c;
    m[11] = 0.0;

    m[12] = 0.0;
    m[13] = 0.0;
    m[14] = 0.0;
    m[15] = 1.0;

    return new Transform(m, Matrix4x4.Transpose(m));
  }

  static Transform LookAt(Point pos, Point look, Vector up) {
    Matrix4x4 m = new Matrix4x4();
    // Initialize fourth column of viewing matrix
    m[3] = pos.x;
    m[7] = pos.y;
    m[11] = pos.z;
    m[15] = 1.0;

    // Initialize first three columns of viewing matrix
    Vector dir = Vector.Normalize(look - pos);
    Vector left = Vector.Normalize(Vector.Cross(Vector.Normalize(up), dir));
    Vector newUp = Vector.Cross(dir, left);
    m[0] = left.x;
    m[4] = left.y;
    m[8] = left.z;
    m[12] = 0.0;
    m[1] = newUp.x;
    m[5] = newUp.y;
    m[9] = newUp.z;
    m[13] = 0.0;
    m[2] = dir.x;
    m[6] = dir.y;
    m[10] = dir.z;
    m[14] = 0.0;

    return new Transform(Matrix4x4.Inverse(m), m);
  }

  static Transform Orthographic(double znear, double zfar) {
    return Scale(1.0, 1.0, 1.0 / (zfar - znear)) *
           Translate(new Vector(0.0, 0.0, -znear));
  }

  static Transform Perspective(double fov, double znear, double zfar) {
    // Perform projective divide
    Matrix4x4 persp = new Matrix4x4.values(
        1.0, 0.0,           0.0,              0.0,
        0.0, 1.0,           0.0,              0.0,
        0.0, 0.0, zfar / (zfar - znear), -zfar * znear / (zfar - znear),
        0.0, 0.0,           1.0,              0.0);

    // Scale to canonical viewing volume
    double invTanAng = 1.0 / Math.tan(Radians(fov) / 2.0);
    return Scale(invTanAng, invTanAng, 1.0) * new Transform(persp);
  }
}
