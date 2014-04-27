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
 * 3D axis aligned bounding box, primarily used to represent the area
 * encompassing one or more objects.
 */
class BBox {
  Point pMin;
  Point pMax;

  BBox([Point p1, Point p2]) {
    if (p1 == null && p2 == null) {
      pMin = new Point(INFINITY, INFINITY, INFINITY);
      pMax = new Point(-INFINITY, -INFINITY, -INFINITY);
    } else if (p1 != null && p2 != null) {
      pMin = new Point(Math.min(p1.x, p2.x),
                       Math.min(p1.y, p2.y),
                       Math.min(p1.z, p2.z));
      pMax = new Point(Math.max(p1.x, p2.x),
                       Math.max(p1.y, p2.y),
                       Math.max(p1.z, p2.z));
    } else {
      pMin = new Point.from(p1);
      pMax = new Point.from(p1);
    }
  }

  BBox.from(BBox other)
      : pMin = new Point.from(other.pMin),
        pMax = new Point.from(other.pMax);

  void reset() {
    pMin.x = INFINITY;
    pMin.y = INFINITY;
    pMin.z = INFINITY;
    pMax.x = -INFINITY;
    pMax.y = -INFINITY;
    pMax.z = -INFINITY;
  }

  void copy(BBox other) {
    pMin.copy(other.pMin);
    pMax.copy(other.pMax);
  }

  double boundingSphere(Point c) {
    c.copy(pMin * 0.5 + pMax * 0.5);
    return inside(c) ? Vector.Distance(c, pMax) : 0.0;
  }

  Point get center => (pMin * 0.5) + (pMax * 0.5);

  Point operator [](int index) => (index == 0) ? pMin : pMax;

  /**
   * Test for a ray intersection with this box. If an intersection was found,
   * true is returned and [hitt0] is set to the closest intersection point
   * (where the ray enters the box) and [hitt1] is set to the farthest
   * intersection point (where the ray exits the box).
   */
  bool intersectP(Ray ray, [List<double> hitt0, List<double> hitt1]) {
    double t0 = ray.minDistance;
    double t1 = ray.maxDistance;
    for (int i = 0; i < 3; ++i) {
      // Update interval for i'th bounding box slab
      double invRayDir = 1.0 / ray.direction[i];
      double tNear = (pMin[i] - ray.origin[i]) * invRayDir;
      double tFar  = (pMax[i] - ray.origin[i]) * invRayDir;

      // Update parametric interval from slab intersection $t$s
      if (tNear > tFar) {
        double t = tNear;
        tNear = tFar;
        tFar = t;
      }

      t0 = tNear > t0 ? tNear : t0;
      t1 = tFar < t1 ? tFar : t1;

      if (t0 > t1) {
        return false;
      }
    }

    if (hitt0 != null) {
      hitt0[0] = t0;
    }

    if (hitt1 != null) {
      hitt1[0] = t1;
    }

    return true;
  }

  bool overlaps(BBox b) {
    bool x = (pMax.x >= b.pMin.x) && (pMin.x <= b.pMax.x);
    bool y = (pMax.y >= b.pMin.y) && (pMin.y <= b.pMax.y);
    bool z = (pMax.z >= b.pMin.z) && (pMin.z <= b.pMax.z);
    return (x && y && z);
  }

  bool inside(Point pt) {
      return (pt.x >= pMin.x && pt.x <= pMax.x &&
              pt.y >= pMin.y && pt.y <= pMax.y &&
              pt.z >= pMin.z && pt.z <= pMax.z);
  }

  BBox setPoint(Point p) {
    pMin.copy(p);
    pMax.copy(p);
    return this;
  }

  BBox unionPoint(Point p) {
    pMin.x = Math.min(pMin.x, p.x);
    pMin.y = Math.min(pMin.y, p.y);
    pMin.z = Math.min(pMin.z, p.z);

    pMax.x = Math.max(pMax.x, p.x);
    pMax.y = Math.max(pMax.y, p.y);
    pMax.z = Math.max(pMax.z, p.z);
    return this;
  }

  BBox union(BBox b2) {
    pMin.x = Math.min(pMin.x, b2.pMin.x);
    pMin.y = Math.min(pMin.y, b2.pMin.y);
    pMin.z = Math.min(pMin.z, b2.pMin.z);

    pMax.x = Math.max(pMax.x, b2.pMax.x);
    pMax.y = Math.max(pMax.y, b2.pMax.y);
    pMax.z = Math.max(pMax.z, b2.pMax.z);
    return this;
  }

  void expand(double delta) {
    pMin.x -= delta;
    pMin.y -= delta;
    pMin.z -= delta;
    pMax.x += delta;
    pMax.y += delta;
    pMax.z += delta;
  }

  double surfaceArea() {
    Vector d = pMax - pMin;
    return 2.0 * (d.x * d.y + d.x * d.z + d.y * d.z);
  }

  double volume() {
    Vector d = pMax - pMin;
    return d.x * d.y * d.z;
  }

  int maximumExtent() {
    Vector diag = pMax - pMin;
    if (diag.x > diag.y && diag.x > diag.z) {
      return 0;
    } else if (diag.y > diag.z) {
      return 1;
    } else {
      return 2;
    }
  }

  Point lerp(double tx, double ty, double tz) {
    return new Point(Lerp(tx, pMin.x, pMax.x),
                     Lerp(ty, pMin.y, pMax.y),
                     Lerp(tz, pMin.z, pMax.z));
  }

  Vector offset(Point p) {
    return new Vector((p.x - pMin.x) / (pMax.x - pMin.x),
                      (p.y - pMin.y) / (pMax.y - pMin.y),
                      (p.z - pMin.z) / (pMax.z - pMin.z));
  }

  static BBox UnionPoint(BBox b, Point p) {
    return new BBox.from(b).unionPoint(p);
  }

  static BBox Union(BBox b, BBox b2) {
    return new BBox.from(b).union(b2);
  }
}
