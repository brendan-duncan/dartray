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

class Disk extends Shape {
  Disk(Transform o2w, Transform w2o, bool ro, this.height,
            this.radius, this.innerRadius, this.phiMax) :
    super(o2w, w2o, ro) {
    phiMax = Radians(phiMax.clamp(0.0, 360.0));
  }

  BBox objectBound() {
    return new BBox(new Point(-radius, -radius, height),
                    new Point( radius,  radius, height));
  }

  static Ray _ray = new Ray();

  bool intersect(Ray r, List<double> tHit, List<double> rayEpsilon,
                 DifferentialGeometry dg) {
    // Transform _Ray_ to object space
    worldToObject.transformRay(r, _ray);

    // Compute plane intersection for disk
    if (_ray.direction.z.abs() < 1.0e-7) {
      return false;
    }

    double thit = (height - _ray.origin.z) / _ray.direction.z;
    if (thit < _ray.minDistance || thit > _ray.maxDistance) {
      return false;
    }

    // See if hit point is inside disk radii and $\phimax$
    Point phit = _ray.pointAt(thit);
    double dist2 = phit.x * phit.x + phit.y * phit.y;
    if (dist2 > radius * radius || dist2 < innerRadius * innerRadius) {
      return false;
    }

    // Test disk $\phi$ value against $\phimax$
    double phi = Math.atan2(phit.y, phit.x);
    if (phi < 0) {
      phi += 2.0 * Math.PI;
    }

    if (phi > phiMax) {
      return false;
    }

    // Find parametric representation of disk hit
    double u = phi / phiMax;
    double oneMinusV = ((Math.sqrt(dist2) - innerRadius) /
                        (radius - innerRadius));

    double invOneMinusV = (oneMinusV > 0.0) ? (1.0 / oneMinusV) : 0.0;
    double v = 1.0 - oneMinusV;
    Vector dpdu = new Vector(-phiMax * phit.y, phiMax * phit.x, 0.0);
    Vector dpdv = new Vector(-phit.x * invOneMinusV, -phit.y * invOneMinusV,
                             0.0);
    dpdu *= phiMax * INV_TWOPI;
    dpdv *= (radius - innerRadius) / radius;
    Normal dndu = new Normal();
    Normal dndv = new Normal();

    // Initialize DifferentialGeometry from parametric information
    Transform o2w = objectToWorld;

    dg.set(o2w.transformPoint(phit),
           o2w.transformVector(dpdu),
           o2w.transformVector(dpdv),
           o2w.transformNormal(dndu),
           o2w.transformNormal(dndv),
           u, v, this);

    // Update tHit for quadric intersection
    tHit[0] = thit;

    // Compute rayEpsilon for quadric intersection
    rayEpsilon[0] = 5.0e-4 * thit;

    return true;
  }

  bool intersectP(Ray r) {
    // Transform Ray to object space
    Ray ray = new Ray();
    worldToObject.transformRay(r, ray);

    // Compute plane intersection for disk
    if (ray.direction.z.abs() < 1.0e-7) {
      return false;
    }

    double thit = (height - ray.origin.z) / ray.direction.z;
    if (thit < ray.minDistance || thit > ray.maxDistance) {
      return false;
    }

    // See if hit point is inside disk radii and phimax
    Point phit = ray.pointAt(thit);
    double dist2 = phit.x * phit.x + phit.y * phit.y;
    if (dist2 > radius * radius || dist2 < innerRadius * innerRadius) {
      return false;
    }

    // Test disk phi value against phimax
    double phi = Math.atan2(phit.y, phit.x);

    if (phi < 0.0) {
      phi += 2.0 * Math.PI;
    }

    if (phi > phiMax) {
      return false;
    }

    return true;
  }

  double area() {
    return phiMax * 0.5 *
           (radius * radius - innerRadius * innerRadius);
  }

  Point sample(double u1, double u2, Normal Ns) {
    List<double> t0 = [0.0];
    List<double> t1 = [0.0];
    ConcentricSampleDisk(u1, u2, t0, t1);
    Point p = new Point(t0[0] * radius, t1[0] * radius, height);
    Ns.copy(objectToWorld.transformNormal(new Normal(0.0, 0.0, 1.0)));
    Ns.normalize();
    if (reverseOrientation) {
      Ns.scale(-1.0);
    }
    return objectToWorld.transformPoint(p);
  }

  double height;
  double radius;
  double innerRadius;
  double phiMax;

  static Disk Create(Transform o2w, Transform w2o,
                     bool reverseOrientation, ParamSet params) {
    double height = params.findOneFloat('height', 0.0);
    double radius = params.findOneFloat('radius', 1.0);
    double inner_radius = params.findOneFloat('innerradius', 0.0);
    double phimax = params.findOneFloat('phimax', 360.0);
    return new Disk(o2w, w2o, reverseOrientation, height, radius,
                         inner_radius, phimax);
  }
}
