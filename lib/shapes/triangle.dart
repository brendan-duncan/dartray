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
part of shapes;

/**
 * An individual triangle in a [TriangleMesh].
 */
class Triangle extends Shape {
  Triangle(Transform o2w, Transform w2o, bool ro, this.mesh, this.index) :
    super(o2w, w2o, ro) {
    index *= 3;
  }

  BBox objectBound() {
    List<Point> tri = mesh.triangle(v(0), v(1), v(2));
    return BBox.UnionPoint(new BBox(worldToObject.transformPoint(tri[0]),
                                    worldToObject.transformPoint(tri[1])),
                           worldToObject.transformPoint(tri[2]));
  }

  BBox worldBound() {
    List<Point> tri = mesh.triangle(v(0), v(1), v(2));
    return BBox.UnionPoint(new BBox(tri[0], tri[1]), tri[2]);
  }

  bool intersect(Ray ray, List<double> tHit, List<double> rayEpsilon,
                 DifferentialGeometry dg) {
    // Get triangle vertices in p1, p2, and p3
    List<Point> tri = mesh.triangle(v(0), v(1), v(2));
    Point p1 = tri[0];
    Point p2 = tri[1];
    Point p3 = tri[2];

    double e1x = p2.data[0] - p1.data[0];
    double e1y = p2.data[1] - p1.data[1];
    double e1z = p2.data[2] - p1.data[2];

    double e2x = p3.data[0] - p1.data[0];
    double e2y = p3.data[1] - p1.data[1];
    double e2z = p3.data[2] - p1.data[2];

    double s1x = (ray.direction.y * e2z) - (ray.direction.z * e2y);
    double s1y = (ray.direction.z * e2x) - (ray.direction.x * e2z);
    double s1z = (ray.direction.x * e2y) - (ray.direction.y * e2x);

    double divisor = (s1x * e1x) + (s1y * e1y) + (s1z * e1z);

    if (divisor == 0.0) {
      return false;
    }

    double invDivisor = 1.0 / divisor;

    // Compute first barycentric coordinate
    double sx = ray.origin.x - p1.x;
    double sy = ray.origin.y - p1.y;
    double sz = ray.origin.z - p1.z;

    double b1 = (sx * s1x + sy * s1y + sz * s1z) * invDivisor;
    if (b1 < 0.0 || b1 > 1.0) {
      return false;
    }

    // Compute second barycentric coordinate
    double s2x = (sy * e1z) - (sz * e1y);
    double s2y = (sz * e1x) - (sx * e1z);
    double s2z = (sx * e1y) - (sy * e1x);

    double b2 = ((ray.direction.x * s2x) + (ray.direction.y * s2y) +
                 (ray.direction.z * s2z)) * invDivisor;

    if (b2 < 0.0 || b1 + b2 > 1.0) {
      return false;
    }

    // Compute t to intersection point
    double t = (e2x * s2x + e2y * s2y + e2z * s2z) * invDivisor;
    if (t < ray.minDistance || t > ray.maxDistance) {
      return false;
    }

    // Compute triangle partial derivatives
    Vector dpdu;
    Vector dpdv;

    getUVs(_uvs);

    // Compute deltas for triangle partial derivatives
    double du1 = _uvs[0] - _uvs[4];
    double du2 = _uvs[2] - _uvs[4];
    double dv1 = _uvs[1] - _uvs[5];
    double dv2 = _uvs[3] - _uvs[5];
    Vector dp1 = p1 - p3;
    Vector dp2 = p2 - p3;
    double determinant = du1 * dv2 - dv1 * du2;

    if (determinant == 0.0) {
      dpdu = new Vector();
      dpdv = new Vector();

      // Cross(e2, e3)
      double e3x = (e2y * e1z) - (e2z * e1y);
      double e3y = (e2z * e1x) - (e2x * e1z);
      double e3z = (e2x * e1y) - (e2y * e1x);
      double len = Math.sqrt(e3x * e3x + e3y * e3y + e3z * e3z);

      // Handle zero determinant for triangle partial derivative matrix
      Vector.CoordinateSystem(new Vector(e3x / len, e3y / len, e3z / len),
                              dpdu, dpdv);
    } else {
      double invdet = 1.0 / determinant;
      dpdu = (dp1 * dv2 - dp2 * dv1) * invdet;
      dpdv = (dp1 * -du2 + dp2 * du1) * invdet;
    }

    // Interpolate (u,v) triangle parametric coordinates
    double b0 = 1.0 - b1 - b2;
    double tu = b0 * _uvs[0] + b1 * _uvs[2] + b2 * _uvs[4];
    double tv = b0 * _uvs[1] + b1 * _uvs[3] + b2 * _uvs[5];

    // Test intersection against alpha texture, if present
    if (ray.depth != -1) {
      if (mesh.alphaTexture != null) {
          DifferentialGeometry dgLocal =
              new DifferentialGeometry().set(ray.pointAt(t), dpdu, dpdv,
                                             Normal.ZERO, Normal.ZERO,
                                             tu, tv, this);

        if (mesh.alphaTexture.evaluate(dgLocal) == 0.0) {
          return false;
        }
      }
    }

    // Fill in DifferentialGeometry from triangle hit
    dg.set(ray.pointAt(t), dpdu, dpdv, Normal.ZERO, Normal.ZERO, tu, tv, this);

    tHit[0] = t;
    rayEpsilon[0] = 1.0e-3 * t;

    return true;
  }

  bool intersectP(Ray ray) {
    // Get triangle vertices in p1, p2, and p3
    List<Point> tri = mesh.triangle(v(0), v(1), v(2));
    Vector e1 = tri[1] - tri[0];
    Vector e2 = tri[2] - tri[0];
    Vector s1 = Vector.Cross(ray.direction, e2);
    double divisor = Vector.Dot(s1, e1);

    if (divisor == 0.0) {
      return false;
    }

    double invDivisor = 1.0 / divisor;

    // Compute first barycentric coordinate
    Vector s = ray.origin - tri[0];
    double b1 = Vector.Dot(s, s1) * invDivisor;
    if (b1 < 0.0 || b1 > 1.0) {
      return false;
    }

    // Compute second barycentric coordinate
    Vector s2 = Vector.Cross(s, e1);
    double b2 = Vector.Dot(ray.direction, s2) * invDivisor;
    if (b2 < 0.0 || b1 + b2 > 1.0) {
      return false;
    }

    // Compute t to intersection point
    double t = Vector.Dot(e2, s2) * invDivisor;
    if (t < ray.minDistance || t > ray.maxDistance) {
      return false;
    }

    // Test shadow ray intersection against alpha texture, if present
    if (ray.depth != -1 && mesh.alphaTexture != null) {
      // Compute triangle partial derivatives
      Vector dpdu;
      Vector dpdv;
      getUVs(_uvs);

      // Compute deltas for triangle partial derivatives
      double du1 = _uvs[0] - _uvs[4];
      double du2 = _uvs[2] - _uvs[4];
      double dv1 = _uvs[1] - _uvs[5];
      double dv2 = _uvs[3] - _uvs[5];
      Vector dp1 = tri[0] - tri[2];
      Vector dp2 = tri[1] - tri[2];
      double determinant = du1 * dv2 - dv1 * du2;

      if (determinant == 0.0) {
        dpdu = new Vector();
        dpdv = new Vector();
        // Handle zero determinant for triangle partial derivative matrix
        Vector.CoordinateSystem(Vector.Normalize(Vector.Cross(e2, e1)),
                                dpdu, dpdv);
      } else {
        double invdet = 1.0 / determinant;
        dpdu = (dp1 * dv2 - dp2 * dv1) * invdet;
        dpdv = (dp1 * -du2 + dp2 * du1) * invdet;
      }

      // Interpolate $(u,v)$ triangle parametric coordinates
      double b0 = 1.0 - b1 - b2;
      double tu = b0 * _uvs[0] + b1 * _uvs[2] + b2 * _uvs[4];
      double tv = b0 * _uvs[1] + b1 * _uvs[3] + b2 * _uvs[5];

      DifferentialGeometry dgLocal =
          new DifferentialGeometry().set(ray.pointAt(t), dpdu, dpdv,
                                         Normal.ZERO, Normal.ZERO,
                                         tu, tv, this);

      if (mesh.alphaTexture.evaluate(dgLocal) == 0.0) {
        return false;
      }
    }

    return true;
  }

  int v(int i) => mesh.vertexIndex[index + i];

  // cache storage for uvs
  static List<double> _uvs = new List<double>(6);

  void getUVs(List<double> uv) {
    if (mesh.uvs != null) {
      uv[0] = mesh.uvs[2 * v(0)];
      uv[1] = mesh.uvs[2 * v(0) + 1];
      uv[2] = mesh.uvs[2 * v(1)];
      uv[3] = mesh.uvs[2 * v(1) + 1];
      uv[4] = mesh.uvs[2 * v(2)];
      uv[5] = mesh.uvs[2 * v(2) + 1];
    } else {
      uv[0] = 0.0;
      uv[1] = 0.0;
      uv[2] = 1.0;
      uv[3] = 0.0;
      uv[4] = 1.0;
      uv[5] = 1.0;
    }
  }

  double area() {
    // Get triangle vertices in _p1_, _p2_, and _p3_
    List<Point> tri = mesh.triangle(v(0), v(1), v(2));
    return 0.5 * Vector.Cross(tri[1] - tri[0], tri[2] - tri[0]).length();
  }

  void getShadingGeometry(Transform obj2world, DifferentialGeometry dg,
                          DifferentialGeometry dgShading) {
    if (mesh.n == null && mesh.s == null) {
      dgShading.copy(dg);
      return;
    }

    final List<double> _by = [0.0];
    final List<double> _bz = [0.0];
    final List<double> _uv = new List<double>(6);

    // Initialize _A_ and _C_ matrices for barycentrics
    getUVs(_uv);

    List<double> A = [_uv[2] - _uv[0],
                      _uv[4] - _uv[0],
                      _uv[3] - _uv[1],
                      _uv[5] - _uv[1]];

    List<double> C = [dg.u - _uv[0], dg.v - _uv[1]];

    double bx;
    if (!SolveLinearSystem2x2(A, C, _by, _bz)) {
      // Handle degenerate parametric mapping
      bx = _by[0] = _bz[0] = 1.0 / 3.0;
    } else {
      bx = 1.0 - _by[0] - _bz[0];
    }

    // Use _n_ and _s_ to compute shading tangents for triangle, _ss_ and _ts_
    Normal ns;
    Vector ss, ts;
    if (mesh.n != null) {
      ns = Vector.Normalize(obj2world.transformNormal(mesh.n[v(0)] * bx +
                                                      mesh.n[v(1)] * _by[0] +
                                                      mesh.n[v(2)] * _bz[0]));
    } else {
      ns = new Normal.from(dg.nn);
    }

    if (mesh.s != null) {
      ss = Vector.Normalize(obj2world.transformVector(mesh.s[v(0)] * bx +
                                                      mesh.s[v(1)] * _by[0] +
                                                      mesh.s[v(2)] * _bz[0]));
    } else {
      ss = Vector.Normalize(dg.dpdu);
    }

    ts = Vector.Cross(ss, ns);
    if (ts.lengthSquared() > 0.0) {
      ts = Vector.Normalize(ts);
      ss = Vector.Cross(ts, ns);
    } else {
      Vector.CoordinateSystem(ns, ss, ts);
    }

    Normal dndu;
    Normal dndv;

    // Compute dndu and dndv for triangle shading geometry
    if (mesh.n != null) {
      // Compute deltas for triangle partial derivatives of normal
      double du1 = _uv[0] - _uv[4];
      double du2 = _uv[2] - _uv[4];
      double dv1 = _uv[1] - _uv[5];
      double dv2 = _uv[3] - _uv[5];
      Normal dn1 = mesh.n[v(0)] - mesh.n[v(2)];
      Normal dn2 = mesh.n[v(1)] - mesh.n[v(2)];
      double determinant = du1 * dv2 - dv1 * du2;
      if (determinant == 0.0) {
        dndu = new Normal();
        dndv = new Normal();
      } else {
        double invdet = 1.0 / determinant;
        dndu = (dn1 * dv2 - dn2 * dv1) * invdet;
        dndv = (dn1 * -du2 + dn2 * du1) * invdet;
      }
    } else {
      dndu = new Normal();
      dndv = new Normal();
    }

    dgShading.set(dg.p, ss, ts,
                 objectToWorld.transformNormal(dndu),
                 objectToWorld.transformNormal(dndv),
                 dg.u, dg.v, dg.shape);

    dgShading.dudx = dg.dudx;
    dgShading.dvdx = dg.dvdx;
    dgShading.dudy = dg.dudy;
    dgShading.dvdy = dg.dvdy;
    dgShading.dpdx = dg.dpdx;
    dgShading.dpdy = dg.dpdy;
  }

  Point sample(double u1, double u2, Normal Ns) {
    List<double> b1 = [0.0];
    List<double> b2 = [0.0];
    UniformSampleTriangle(u1, u2, b1, b2);

    // Get triangle vertices in _p1_, _p2_, and _p3_
    List<Point> tri = mesh.triangle(v(0), v(1), v(2));
    Point p = tri[0] * b1[0] + tri[1] * b2[0] + tri[2] * (1.0 - b1[0] - b2[0]);
    Vector n = Vector.Cross(tri[1] - tri[0], tri[2] - tri[0]);
    Ns.copy(Vector.Normalize(n));
    if (reverseOrientation) {
      Ns.x *= -1.0;
      Ns.y *= -1.0;
      Ns.z *= -1.0;
    }

    return p;
  }

  TriangleMesh mesh;
  int index;
}
