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

class Hyperboloid extends Shape {
  Hyperboloid(Transform o2w, Transform w2o, bool ro,
                   this.p1, this.p2, double tm) :
    super(o2w, w2o, ro) {
    double radius1 = Math.sqrt(p1.x * p1.x + p1.y * p1.y);
    double radius2 = Math.sqrt(p2.x * p2.x + p2.y * p2.y);
    rmax = Math.max(radius1, radius2);
    zmin = Math.min(p1.z, p2.z);
    zmax = Math.max(p1.z, p2.z);
    phiMax = Radians(tm.clamp(0.0, 360.0));

    // Compute implicit function coefficients for hyperboloid
    if (p2.z == 0.0) {
      Point t = p1;
      p1 = p2;
      p2 = t;
    }

     Point pp = p1;
     double xy1;
     double xy2;
     do {
         pp += (p2 - p1) * 2.0;
         xy1 = pp.x * pp.x + pp.y * pp.y;
         xy2 = p2.x * p2.x + p2.y * p2.y;
         a = (1.0 / xy1 - (pp.z * pp.z) / (xy1 * p2.z * p2.z)) /
             (1 - (xy2 * pp.z * pp.z) / (xy1 * p2.z * p2.z));
         c = (a * xy2 - 1) / (p2.z * p2.z);
     } while (a.isInfinite || a.isNaN);
  }

  BBox objectBound() {
    Point p1 = new Point(-rmax, -rmax, zmin);
    Point p2 = new Point(rmax,  rmax, zmax);
    return new BBox(p1, p2);
  }

  bool intersect(Ray r, List<double> tHit, List<double> rayEpsilon,
                 DifferentialGeometry dg) {
    // Transform _Ray_ to object space
    Ray ray = new Ray();
    worldToObject.transformRay(r, ray);

    // Compute quadratic hyperboloid coefficients
    double A = a * ray.direction.x * ray.direction.x +
               a * ray.direction.y * ray.direction.y -
               c * ray.direction.z * ray.direction.z;
    double B = 2.0 * (a * ray.direction.x * ray.origin.x +
                      a * ray.direction.y * ray.origin.y -
                      c * ray.direction.z * ray.origin.z);
    double C = a * ray.origin.x * ray.origin.x +
               a * ray.origin.y * ray.origin.y -
               c * ray.origin.z * ray.origin.z - 1;

    // Solve quadratic equation for _t_ values
    List<double> t0 = [0.0];
    List<double> t1 = [0.0];
    if (!Quadratic(A, B, C, t0, t1)) {
      return false;
    }

    // Compute intersection distance along ray
    if (t0[0] > ray.maxDistance || t1[0] < ray.minDistance) {
      return false;
    }

    double thit = t0[0];
    if (t0[0] < ray.minDistance) {
      thit = t1[0];
      if (thit > ray.maxDistance) {
        return false;
      }
    }

    // Compute hyperboloid inverse mapping
    Point phit = ray.pointAt(thit);
    double v = (phit.z - p1.z) / (p2.z - p1.z);
    Point pr = p1 * (1.0 - v) + p2 * v;
    double phi = Math.atan2(pr.x * phit.y - phit.x * pr.y,
                            phit.x * pr.x + phit.y * pr.y);
    if (phi < 0) {
      phi += 2.0 * Math.PI;
    }

    // Test hyperboloid intersection against clipping parameters
    if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
      if (thit == t1[0]) {
        return false;
      }
      thit = t1[0];
      if (t1[0] > ray.maxDistance) {
        return false;
      }
      // Compute hyperboloid inverse mapping
      phit = ray.pointAt(thit);
      v = (phit.z - p1.z) / (p2.z - p1.z);
      Point pr = p1 * (1.0 - v) + p2 * v;
      phi = Math.atan2(pr.x * phit.y - phit.x * pr.y,
                       phit.x * pr.x + phit.y * pr.y);
      if (phi < 0.0) {
        phi += 2.0 * Math.PI;
      }

      if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
        return false;
      }
    }

    // Compute parametric representation of hyperboloid hit
    double u = phi / phiMax;

    // Compute hyperboloid $\dpdu$ and $\dpdv$
    double cosphi = Math.cos(phi);
    double sinphi = Math.sin(phi);
    Vector dpdu = new Vector(-phiMax * phit.y, phiMax * phit.x, 0.0);
    Vector dpdv = new Vector((p2.x - p1.x) * cosphi - (p2.y - p1.y) * sinphi,
                             (p2.x - p1.x) * sinphi + (p2.y - p1.y) * cosphi,
                             p2.z - p1.z);

    // Compute hyperboloid $\dndu$ and $\dndv$
    Vector d2Pduu = new Vector(phit.x, phit.y, 0.0) * (-phiMax * phiMax);
    Vector d2Pduv = new Vector(-dpdv.y, dpdv.x, 0.0) * phiMax;
    Vector d2Pdvv = new Vector();

    // Compute coefficients for fundamental forms
    double E = Vector.Dot(dpdu, dpdu);
    double F = Vector.Dot(dpdu, dpdv);
    double G = Vector.Dot(dpdv, dpdv);
    Vector N = Vector.Normalize(Vector.Cross(dpdu, dpdv));
    double e = Vector.Dot(N, d2Pduu);
    double f = Vector.Dot(N, d2Pduv);
    double g = Vector.Dot(N, d2Pdvv);

    // Compute $\dndu$ and $\dndv$ from fundamental form coefficients
    double invEGF2 = 1.0 / (E * G - F * F);
    Normal dndu = new Normal.from(dpdu * ((f * F - e * G) * invEGF2) +
                                  dpdv * ((e * F - f * E) * invEGF2));
    Normal dndv = new Normal.from(dpdu * ((g * F - f * G) * invEGF2) +
                                  dpdv * ((f * F - g * E) * invEGF2));

    // Initialize _DifferentialGeometry_ from parametric information
    Transform o2w = objectToWorld;
    dg.set(o2w.transformPoint(phit),
           o2w.transformVector(dpdu),
           o2w.transformVector(dpdv),
           o2w.transformNormal(dndu),
           o2w.transformNormal(dndv),
           u, v, this);

    // Update _tHit_ for quadric intersection
    tHit[0] = thit;

    // Compute _rayEpsilon_ for quadric intersection
    rayEpsilon[0] = 5.0e-4 * thit;

    return true;
  }

  bool intersectP(Ray r) {
    // Transform _Ray_ to object space
    Ray ray = new Ray();
    worldToObject.transformRay(r, ray);

    // Compute quadratic hyperboloid coefficients
    double A = a * ray.direction.x * ray.direction.x +
               a * ray.direction.y * ray.direction.y -
               c * ray.direction.z * ray.direction.z;
    double B = 2.0 * (a * ray.direction.x * ray.origin.x +
                      a * ray.direction.y * ray.origin.y -
                      c * ray.direction.z * ray.origin.z);
    double C = a * ray.origin.x * ray.origin.x +
               a * ray.origin.y * ray.origin.y -
               c * ray.origin.z * ray.origin.z - 1;

    // Solve quadratic equation for _t_ values
    List<double> t0 = [0.0];
    List<double> t1 = [0.0];
    if (!Quadratic(A, B, C, t0, t1)) {
      return false;
    }

    // Compute intersection distance along ray
    if (t0[0] > ray.maxDistance || t1[0] < ray.minDistance) {
      return false;
    }

    double thit = t0[0];
    if (t0[0] < ray.minDistance) {
      thit = t1[0];
      if (thit > ray.maxDistance) {
        return false;
      }
    }

    // Compute hyperboloid inverse mapping
    Point phit = ray.pointAt(thit);
    double v = (phit.z - p1.z) / (p2.z - p1.z);
    Point pr = p1 * (1.0 - v) + p2 * v;
    double phi = Math.atan2(pr.x * phit.y - phit.x * pr.y,
                            phit.x * pr.x + phit.y * pr.y);
    if (phi < 0) {
      phi += 2.0 * Math.PI;
    }

    // Test hyperboloid intersection against clipping parameters
    if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
      if (thit == t1[0]) {
        return false;
      }
      thit = t1[0];
      if (t1[0] > ray.maxDistance) {
        return false;
      }
      // Compute hyperboloid inverse mapping
      phit = ray.pointAt(thit);
      v = (phit.z - p1.z) / (p2.z - p1.z);
      Point pr = p1 * (1.0 - v) + p2 * v;
      phi = Math.atan2(pr.x * phit.y - phit.x * pr.y,
                       phit.x * pr.x + phit.y * pr.y);
      if (phi < 0.0) {
        phi += 2.0 * Math.PI;
      }

      if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
        return false;
      }
    }

    return true;
  }

  double area() {
    SQR(a) => a * a;
    QUAD(a) => a * a * a * a;
    return phiMax / 6.0 *
         (2.0 * QUAD(p1.x) - 2.0 * p1.x * p1.x * p1.x * p2.x +
          2.0 * QUAD(p2.x) + 2.0 * (p1.y * p1.y + p1.y * p2.y + p2.y * p2.y) *
          (SQR(p1.y - p2.y) + SQR(p1.z - p2.z)) +
          p2.x * p2.x * (5.0 * p1.y * p1.y + 2.0 * p1.y * p2.y -
          4.0 * p2.y * p2.y + 2.0 * SQR(p1.z - p2.z)) +
          p1.x * p1.x*(-4.0 * p1.y * p1.y + 2.0 * p1.y * p2.y +
          5.0 * p2.y * p2.y + 2.0 * SQR(p1.z - p2.z)) -
          2.0 * p1.x * p2.x * (p2.x * p2.x - p1.y * p1.y +
          5.0 * p1.y * p2.y - p2.y * p2.y - p1.z * p1.z +
          2.0 * p1.z * p2.z - p2.z * p2.z));
  }

  Point p1;
  Point p2;
  double zmin;
  double zmax;
  double phiMax;
  double rmax;
  double a;
  double c;

  static Hyperboloid Create(Transform o2w, Transform w2o,
                            bool reverseOrientation, ParamSet params) {
    Point p1 = params.findOnePoint('p1', new Point(0.0, 0.0, 0.0));
    Point p2 = params.findOnePoint('p2', new Point(1.0, 1.0, 1.0));
    double phimax = params.findOneFloat('phimax', 360.0);
    return new Hyperboloid(o2w, w2o, reverseOrientation, p1, p2, phimax);
  }
}
