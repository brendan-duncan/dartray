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

class Anisotropic extends MicrofacetDistribution {
  double ex;
  double ey;

  Anisotropic(this.ex, this.ey) {
    if (ex > 10000.0 || ex.isNaN) {
      ex = 10000.0;
    }
    if (ey > 10000.0 || ey.isNaN) {
      ey = 10000.0;
    }
  }

  double d(Vector wh) {
    double costhetah = wh.z.abs();
    double d = 1.0 - costhetah * costhetah;
    if (d == 0.0) {
      return 0.0;
    }

    double e = (ex * wh.x * wh.x + ey * wh.y * wh.y) / d;

    return Math.sqrt((ex + 2.0) * (ey + 2.0)) *
           INV_TWOPI * Math.pow(costhetah, e);
  }

  double sample_f(Vector wo, Vector wi, double u1, double u2) {
    // Sample from first quadrant and remap to hemisphere to sample wh
    List<double> phi_cosTheta = [0.0, 0.0];
    if (u1 < 0.25) {
      sampleFirstQuadrant(4.0 * u1, u2, phi_cosTheta);
    } else if (u1 < 0.5) {
      u1 = 4.0 * (0.5 - u1);
      sampleFirstQuadrant(u1, u2, phi_cosTheta);
      phi_cosTheta[0] = Math.PI - phi_cosTheta[0];
    } else if (u1 < 0.75) {
      u1 = 4.0 * (u1 - 0.5);
      sampleFirstQuadrant(u1, u2, phi_cosTheta);
      phi_cosTheta[0] += Math.PI;
    } else {
      u1 = 4.0 * (1.0 - u1);
      sampleFirstQuadrant(u1, u2, phi_cosTheta);
      phi_cosTheta[0] = 2.0 * Math.PI - phi_cosTheta[0];
    }

    double phi = phi_cosTheta[0];
    double cosTheta = phi_cosTheta[1];

    double sintheta = Math.sqrt(Math.max(0.0, 1.0 - cosTheta * cosTheta));
    Vector wh = Vector.SphericalDirection(sintheta, cosTheta, phi);
    if (!Vector.SameHemisphere(wo, wh)) {
      wh = -wh;
    }

    // Compute incident direction by reflecting about wh
    wi.copy(-wo + wh * 2.0 * Vector.Dot(wo, wh));

    // Compute PDF for $\wi$ from anisotropic distribution
    double costhetah = Vector.AbsCosTheta(wh);
    double ds = 1.0 - costhetah * costhetah;
    double anisotropic_pdf = 0.0;
    if (ds > 0.0 && Vector.Dot(wo, wh) > 0.0) {
       double e = (ex * wh.x * wh.x + ey * wh.y * wh.y) / ds;
       double d = Math.sqrt((ex + 1.0) * (ey + 1.0)) * INV_TWOPI *
                 Math.pow(costhetah, e);
       anisotropic_pdf = d / (4.0 * Vector.Dot(wo, wh));
    }

    return anisotropic_pdf;
  }

  double pdf(Vector wo, Vector wi) {
    Vector wh = Vector.Normalize(wo + wi);
    // Compute PDF for [wi] from anisotropic distribution
    double costhetah = Vector.AbsCosTheta(wh);
    double ds = 1.0 - costhetah * costhetah;
    double anisotropic_pdf = 0.0;
    if (ds > 0.0 && Vector.Dot(wo, wh) > 0.0) {
      double e = (ex * wh.x * wh.x + ey * wh.y * wh.y) / ds;
      double d = Math.sqrt((ex + 1.0) * (ey + 1.0)) * INV_TWOPI *
                 Math.pow(costhetah, e);
      anisotropic_pdf = d / (4.0 * Vector.Dot(wo, wh));
    }
    return anisotropic_pdf;
  }

  void sampleFirstQuadrant(double u1, double u2, List<double> phi_costheta) {
    if (ex == ey) {
      phi_costheta[0] = Math.PI * u1 * 0.5;
    } else {
      phi_costheta[0] = Math.atan(Math.sqrt((ex + 1.0) / (ey + 1.0)) *
                        Math.tan(Math.PI * u1 * 0.5));
    }

    double cosphi = Math.cos(phi_costheta[0]);
    double sinphi = Math.sin(phi_costheta[0]);

    phi_costheta[1] = Math.pow(u2, 1.0 / (ex * cosphi * cosphi +
                                          ey * sinphi * sinphi + 1.0));
  }
}
