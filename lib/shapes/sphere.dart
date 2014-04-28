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

class Sphere extends Shape {
  Sphere(Transform o2w, Transform w2o, bool ro, this.radius,
         double z0, double z1, double pm) :
    super(o2w, w2o, ro) {
    zmin = Math.min(z0, z1).clamp(-radius, radius);
    zmax = Math.max(z0, z1).clamp(-radius, radius);
    thetaMin = Math.acos((zmin / radius).clamp(-1.0, 1.0));
    thetaMax = Math.acos((zmax / radius).clamp(-1.0, 1.0));
    phiMax = Radians(pm.clamp(0.0, 360.0));
  }

  BBox objectBound() {
    return new BBox(new Point(-radius, -radius, zmin),
                    new Point( radius,  radius, zmax));
  }

  bool intersect(Ray r, List<double> tHit, List<double> rayEpsilon,
                 DifferentialGeometry dg) {
    double phi;
    Point phit = new Point();
    // Transform _Ray_ to object space
    Ray ray = new Ray();
    worldToObject.transformRay(r, ray);

    // Compute quadratic sphere coefficients
    double A = ray.direction.x * ray.direction.x +
               ray.direction.y * ray.direction.y +
               ray.direction.z * ray.direction.z;
    double B = 2 * (ray.direction.x * ray.origin.x +
                    ray.direction.y * ray.origin.y +
                    ray.direction.z * ray.origin.z);
    double C = ray.origin.x * ray.origin.x +
               ray.origin.y * ray.origin.y +
               ray.origin.z * ray.origin.z -
               radius * radius;

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
    if (thit < ray.minDistance) {
      thit = t1[0];
      if (thit > ray.maxDistance) {
        return false;
      }
    }

    // Compute sphere hit position and $\phi$
    phit = ray.pointAt(thit);
    if (phit.x == 0.0 && phit.y == 0.0) {
      phit.x = 1.0e-5 * radius;
    }

    phi = Math.atan2(phit.y, phit.x);
    if (phi < 0.0) {
      phi += 2.0 * Math.PI;
    }

    // Test sphere intersection against clipping parameters
    if ((zmin > -radius && phit.z < zmin) ||
        (zmax <  radius && phit.z > zmax) || phi > phiMax) {
      if (thit == t1[0]) {
        return false;
      }
      if (t1[0] > ray.maxDistance) {
        return false;
      }
      thit = t1[0];

      // Compute sphere hit position and $\phi$
      phit = ray.pointAt(thit);
      if (phit.x == 0.0 && phit.y == 0.0) {
        phit.x = 1.0e-5 * radius;
      }

      phi = Math.atan2(phit.y, phit.x);
      if (phi < 0.0) {
        phi += 2.0 * Math.PI;
      }
      if ((zmin > -radius && phit.z < zmin) ||
          (zmax <  radius && phit.z > zmax) || phi > phiMax) {
        return false;
      }
    }

    // Find parametric representation of sphere hit
    double u = phi / phiMax;
    double theta = Math.acos((phit.z / radius).clamp(-1.0, 1.0));
    double v = (theta - thetaMin) / (thetaMax - thetaMin);

    // Compute sphere $\dpdu$ and $\dpdv$
    double zradius = Math.sqrt(phit.x * phit.x + phit.y * phit.y);
    double invzradius = 1.0 / zradius;
    double cosphi = phit.x * invzradius;
    double sinphi = phit.y * invzradius;
    Vector dpdu = new Vector(-phiMax * phit.y, phiMax * phit.x, 0.0);
    Vector dpdv = new Vector(phit.z * cosphi, phit.z * sinphi,
                             -radius * Math.sin(theta)) * (thetaMax - thetaMin);

    // Compute sphere $\dndu$ and $\dndv$
    Vector d2Pduu = new Vector(phit.x, phit.y, 0.0) * -phiMax * phiMax;
    Vector d2Pduv = new Vector(-sinphi, cosphi, 0.0) *
                    (thetaMax - thetaMin) * phit.z * phiMax;
    Vector d2Pdvv = new Vector(phit.x, phit.y, phit.z) *
                    -(thetaMax - thetaMin) * (thetaMax - thetaMin);

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
    Normal dndu = new Normal.from(dpdu * (f * F - e * G) * invEGF2 +
                                  dpdv * (e * F - f * E) * invEGF2);
    Normal dndv = new Normal.from(dpdu * (g * F - f * G) * invEGF2 +
                                  dpdv * (f * F - g * E) * invEGF2);

    // Initialize _DifferentialGeometry_ from parametric information
    Transform o2w = objectToWorld;
    dg.set(o2w.transformPoint(phit), o2w.transformVector(dpdu),
          o2w.transformVector(dpdv), o2w.transformNormal(dndu),
          o2w.transformNormal(dndv), u, v, this);

    // Update _tHit_ for quadric intersection
    tHit[0] = thit;

    // Compute _rayEpsilon_ for quadric intersection
    rayEpsilon[0] = 5.0e-4 * tHit[0];

    return true;
  }

  bool intersectP(Ray r) {
    double phi;
    Point phit;
    // Transform _Ray_ to object space
    Ray ray = new Ray();
    worldToObject.transformRay(r, ray);

    // Compute quadratic sphere coefficients
    double A = ray.direction.x * ray.direction.x + ray.direction.y * ray.direction.y +
               ray.direction.z * ray.direction.z;
    double B = 2 * (ray.direction.x * ray.origin.x +
                    ray.direction.y * ray.origin.y +
                    ray.direction.z * ray.origin.z);
    double C = ray.origin.x * ray.origin.x + ray.origin.y * ray.origin.y +
               ray.origin.z * ray.origin.z - radius * radius;

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
    if (thit < ray.minDistance) {
      thit = t1[0];
      if (thit > ray.maxDistance) {
        return false;
      }
    }

    // Compute sphere hit position and $\phi$
    phit = ray.pointAt(thit);
    if (phit.x == 0.0 && phit.y == 0.0) {
      phit.x = 1.0e-5 * radius;
    }
    phi = Math.atan2(phit.y, phit.x);
    if (phi < 0.0) {
      phi += 2.0 * Math.PI;
    }

    // Test sphere intersection against clipping parameters
    if ((zmin > -radius && phit.z < zmin) ||
        (zmax <  radius && phit.z > zmax) || phi > phiMax) {
      if (thit == t1) {
        return false;
      }
      if (t1[0] > ray.maxDistance) {
        return false;
      }
      thit = t1[0];
      // Compute sphere hit position and $\phi$
      phit = ray.pointAt(thit);
      if (phit.x == 0.0 && phit.y == 0.0) {
        phit.x = 1.0e-5 * radius;
      }
      phi = Math.atan2(phit.y, phit.x);
      if (phi < 0.0) {
        phi += 2.0 * Math.PI;
      }
      if ((zmin > -radius && phit.z < zmin) ||
          (zmax <  radius && phit.z > zmax) || phi > phiMax) {
        return false;
      }
    }

    return true;
  }

  double area() {
    return phiMax * radius * (zmax - zmin);
  }

  Point sample(double u1, double u2, Normal ns) {
    Point p = new Point() + UniformSampleSphere(u1, u2) * radius;
    ns.copy(Vector.Normalize(objectToWorld.transformNormal(new Normal(p.x, p.y, p.z))));
    if (reverseOrientation) {
      ns.x = -ns.x;
      ns.y = -ns.y;
      ns.z = -ns.z;
    }
    return objectToWorld.transformPoint(p);
  }

  Point sample2(Point p, double u1, double u2, Normal ns) {
    // Compute coordinate system for sphere sampling
    Point Pcenter = objectToWorld.transformPoint(new Point());
    Vector wc = Vector.Normalize(Pcenter - p);
    Vector wcX = new Vector();
    Vector wcY = new Vector();
    Vector.CoordinateSystem(wc, wcX, wcY);

    // Sample uniformly on sphere if $\pt{}$ is inside it
    if (Vector.DistanceSquared(p, Pcenter) - radius * radius < 1.0e-4) {
      return sample(u1, u2, ns);
    }

    // Sample sphere uniformly inside subtended cone
    double sinThetaMax2 = radius * radius / Vector.DistanceSquared(p, Pcenter);
    double cosThetaMax = Math.sqrt(Math.max(0.0, 1.0 - sinThetaMax2));
    DifferentialGeometry dgSphere = new DifferentialGeometry();
    List<double> thit = [0.0];
    List<double> rayEpsilon = [0.0];
    Point ps = new Point();
    Ray r = new Ray(p, UniformSampleCone2(u1, u2, cosThetaMax, wcX, wcY, wc),
                    1.0e-3);

    if (!intersect(r, thit, rayEpsilon, dgSphere)) {
      thit[0] = Vector.Dot(Pcenter - p, Vector.Normalize(r.direction));
    }

    ps = r.pointAt(thit[0]);
    ns.copy(Vector.Normalize(ps - Pcenter));

    if (reverseOrientation) {
      ns.x = -ns.x;
      ns.y = -ns.y;
      ns.z = -ns.z;
    }

    return ps;
  }

  double pdf2(Point p, Vector wi) {
    Point Pcenter = objectToWorld.transformPoint(new Point());
    // Return uniform weight if point inside sphere
    if (Vector.DistanceSquared(p, Pcenter) - radius * radius < 1.0e-4) {
      return super.pdf2(p, wi);
    }

    // Compute general sphere weight
    double sinThetaMax2 = radius * radius / Vector.DistanceSquared(p, Pcenter);
    double cosThetaMax = Math.sqrt(Math.max(0.0, 1.0 - sinThetaMax2));
    return UniformConePdf(cosThetaMax);
  }

  double radius;
  double phiMax;
  double zmin, zmax;
  double thetaMin, thetaMax;

  static Sphere Create(Transform o2w, Transform w2o,
                       bool reverseOrientation, ParamSet params) {
    double radius = params.findOneFloat('radius', 1.0);
    double zmin = params.findOneFloat('zmin', -radius);
    double zmax = params.findOneFloat('zmax', radius);
    double phimax = params.findOneFloat('phimax', 360.0);
    return new Sphere(o2w, w2o, reverseOrientation, radius,
                      zmin, zmax, phimax);
  }
}
