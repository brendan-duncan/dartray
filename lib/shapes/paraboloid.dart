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

class Paraboloid extends Shape {
  Paraboloid(Transform o2w, Transform w2o, bool ro, this.radius,
             double z0, double z1, double tm) :
    super(o2w, w2o, ro) {
    zmin = Math.min(z0, z1);
    zmax = Math.max(z0, z1);
    phiMax = Radians(tm.clamp(0.0, 360.0));
  }

  BBox objectBound() {
    Point p1 = new Point(-radius, -radius, zmin);
    Point p2 = new Point(radius,  radius, zmax);
    return new BBox(p1, p2);
  }

  bool intersect(Ray r, List<double> tHit, List<double> rayEpsilon,
                 DifferentialGeometry dg) {
    // Transform _Ray_ to object space
    Ray ray = new Ray();
    worldToObject.transformRay(r, ray);

    // Compute quadratic paraboloid coefficients
    double k = zmax / (radius * radius);
    double A = k * (ray.direction.x * ray.direction.x +
                    ray.direction.y * ray.direction.y);
    double B = 2 * k * (ray.direction.x * ray.origin.x +
                        ray.direction.y * ray.origin.y) - ray.direction.z;
    double C = k * (ray.origin.x * ray.origin.x +
                    ray.origin.y * ray.origin.y) - ray.origin.z;

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

    // Compute paraboloid inverse mapping
    Point phit = ray.pointAt(thit);
    double phi = Math.atan2(phit.y, phit.x);
    if (phi < 0.0) {
      phi += 2.0 * Math.PI;
    }

    // Test paraboloid intersection against clipping parameters
    if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
      if (thit == t1[0]) {
        return false;
      }

      thit = t1[0];
      if (t1[0] > ray.maxDistance) {
        return false;
      }

      // Compute paraboloid inverse mapping
      phit = ray.pointAt(thit);
      phi = Math.atan2(phit.y, phit.x);
      if (phi < 0.0) {
        phi += 2.0 * Math.PI;
      }

      if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
        return false;
      }
    }

    // Find parametric representation of paraboloid hit
    double u = phi / phiMax;
    double v = (phit.z - zmin) / (zmax - zmin);

    // Compute parabaloid $\dpdu$ and $\dpdv$
    Vector dpdu = new Vector(-phiMax * phit.y, phiMax * phit.x, 0.0);
    Vector dpdv = new Vector(phit.x / (2.0 * phit.z),
                             phit.y / (2.0 * phit.z), 1.0) * (zmax - zmin);

    // Compute parabaloid $\dndu$ and $\dndv$
    Vector d2Pduu = new Vector(phit.x, phit.y, 0.0) * (-phiMax * phiMax);
    Vector d2Pduv = new Vector(-phit.y / (2.0 * phit.z),
                               phit.x / (2.0 * phit.z),
                               0.0) * (zmax - zmin) * phiMax;
    Vector d2Pdvv = new Vector(phit.x / (4.0 * phit.z * phit.z),
                               phit.y / (4.0 * phit.z * phit.z),
                               0.0) * (-(zmax - zmin) * (zmax - zmin));

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

    // Compute quadratic paraboloid coefficients
    double k = zmax / (radius * radius);
    double A = k * (ray.direction.x * ray.direction.x +
                    ray.direction.y * ray.direction.y);
    double B = 2 * k * (ray.direction.x * ray.origin.x +
                        ray.direction.y * ray.origin.y) - ray.direction.z;
    double C = k * (ray.origin.x * ray.origin.x +
                    ray.origin.y * ray.origin.y) - ray.origin.z;

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

    // Compute paraboloid inverse mapping
    Point phit = ray.pointAt(thit);
    double phi = Math.atan2(phit.y, phit.x);
    if (phi < 0.0) {
      phi += 2.0 * Math.PI;
    }

    // Test paraboloid intersection against clipping parameters
    if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
      if (thit == t1[0]) {
        return false;
      }

      thit = t1[0];
      if (t1[0] > ray.maxDistance) {
        return false;
      }

      // Compute paraboloid inverse mapping
      phit = ray.pointAt(thit);
      phi = Math.atan2(phit.y, phit.x);
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
    return phiMax / 12.0 * (Math.pow(1 + 4 * zmin, 1.5) -
           Math.pow(1 + 4 * zmax, 1.5));
  }

  double radius;
  double zmin;
  double zmax;
  double phiMax;

  static Paraboloid Create(Transform o2w, Transform w2o,
                                bool reverseOrientation, ParamSet params) {
      double radius = params.findOneFloat("radius", 1.0);
      double zmin = params.findOneFloat("zmin", 0.0);
      double zmax = params.findOneFloat("zmax", 1.0);
      double phimax = params.findOneFloat("phimax", 360.0);
      return new Paraboloid(o2w, w2o, reverseOrientation,
                                 radius, zmin, zmax, phimax);
  }
}
