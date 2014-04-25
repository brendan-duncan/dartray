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

class Cone extends Shape {
  Cone(Transform o2w, Transform w2o, bool ro,
       this.height, this.radius, double tm) :
    super(o2w, w2o, ro) {
    phiMax = Radians(tm.clamp(0.0, 360.0));
  }

  BBox objectBound() {
    Point p1 = new Point(-radius, -radius, 0.0);
    Point p2 = new Point(radius,  radius, height);
    return new BBox(p1, p2);
  }

  bool intersect(Ray r, List<double> tHit, List<double> rayEpsilon,
                 DifferentialGeometry dg) {
    // Transform _Ray_ to object space
    Ray ray = new Ray();
    worldToObject.transformRay(r, ray);

    // Compute quadratic cone coefficients
    double k = radius / height;
    k = k * k;
    double A = ray.direction.x * ray.direction.x +
               ray.direction.y * ray.direction.y -
               k * ray.direction.z * ray.direction.z;
    double B = 2.0 * (ray.direction.x * ray.origin.x +
                      ray.direction.y * ray.origin.y -
                      k * ray.direction.z * (ray.origin.z - height));
    double C = ray.origin.x * ray.origin.x + ray.origin.y * ray.origin.y -
               k * (ray.origin.z - height) * (ray.origin.z - height);

    // Solve quadratic equation for _t_ values
    List<double> t0 = [0.0],
                 t1 = [0.0];
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

    // Compute cone inverse mapping
    Point phit = ray.pointAt(thit);
    double phi = Math.atan2(phit.y, phit.x);
    if (phi < 0.0) {
      phi += 2.0 * Math.PI;
    }

    // Test cone intersection against clipping parameters
    if (phit.z < 0 || phit.z > height || phi > phiMax) {
      if (thit == t1[0]) {
        return false;
      }
      thit = t1[0];
      if (t1[0] > ray.maxDistance) {
        return false;
      }

      // Compute cone inverse mapping
      phit = ray.pointAt(thit);
      phi = Math.atan2(phit.y, phit.x);
      if (phi < 0.0) {
        phi += 2.0 * Math.PI;
      }
      if (phit.z < 0 || phit.z > height || phi > phiMax) {
        return false;
      }
    }

    // Find parametric representation of cone hit
    double u = phi / phiMax;
    double v = phit.z / height;

    // Compute cone $\dpdu$ and $\dpdv$
    Vector dpdu = new Vector(-phiMax * phit.y, phiMax * phit.x, 0.0);
    Vector dpdv = new Vector(-phit.x / (1.0 - v),
                             -phit.y / (1.0 - v), height);

    // Compute cone $\dndu$ and $\dndv$
    Vector d2Pduu = new Vector(phit.x, phit.y, 0.0) * (-phiMax * phiMax);
    Vector d2Pduv = new Vector(phit.y, -phit.x, 0.0) * (phiMax / (1.0 - v));
    Vector d2Pdvv = new Vector(0.0, 0.0, 0.0);

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

    // Compute quadratic cone coefficients
    double k = radius / height;
    k = k * k;
    double A = ray.direction.x * ray.direction.x +
               ray.direction.y * ray.direction.y -
               k * ray.direction.z * ray.direction.z;
    double B = 2.0 * (ray.direction.x * ray.origin.x +
                      ray.direction.y * ray.origin.y -
                      k * ray.direction.z * (ray.origin.z - height));
    double C = ray.origin.x * ray.origin.x + ray.origin.y * ray.origin.y -
               k * (ray.origin.z - height) * (ray.origin.z - height);

    // Solve quadratic equation for _t_ values
    List<double> t0 = [0.0],
                 t1 = [0.0];
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

    // Compute cone inverse mapping
    Point phit = ray.pointAt(thit);
    double phi = Math.atan2(phit.y, phit.x);
    if (phi < 0.0) {
      phi += 2.0 * Math.PI;
    }

    // Test cone intersection against clipping parameters
    if (phit.z < 0 || phit.z > height || phi > phiMax) {
      if (thit == t1[0]) {
        return false;
      }
      thit = t1[0];
      if (t1[0] > ray.maxDistance) {
        return false;
      }

      // Compute cone inverse mapping
      phit = ray.pointAt(thit);
      phi = Math.atan2(phit.y, phit.x);
      if (phi < 0.0) {
        phi += 2.0 * Math.PI;
      }
      if (phit.z < 0 || phit.z > height || phi > phiMax) {
        return false;
      }
    }

    return true;
  }

  double area() {
    return radius * Math.sqrt((height * height) + (radius * radius)) *
           phiMax / 2.0;
  }

  double radius;
  double height;
  double phiMax;

  static Cone Create(Transform o2w, Transform w2o,
                     bool reverseOrientation, ParamSet params) {
    double radius = params.findOneFloat('radius', 1.0);
    double height = params.findOneFloat('height', 1.0);
    double phimax = params.findOneFloat('phimax', 360.0);
    return new Cone(o2w, w2o, reverseOrientation, height, radius, phimax);
  }
}
