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

class Blinn extends MicrofacetDistribution {
  Blinn(this.exponent) {
    if (exponent > 10000.0 || exponent.isNaN) {
      exponent = 10000.0;
    }
  }

  double d(Vector wh) {
    double costhetah = Vector.AbsCosTheta(wh);
    return (exponent + 2.0) * INV_TWOPI * Math.pow(costhetah, exponent);
  }

  double sample_f(Vector wo, Vector wi, double u1, double u2) {
    // Compute sampled half-angle vector [wh] for Blinn distribution
    double costheta = Math.pow(u1, 1.0 / (exponent + 1.0));
    double sintheta = Math.sqrt(Math.max(0.0, 1.0 - costheta * costheta));
    double phi = u2 * 2.0 * Math.PI;
    Vector wh = Vector.SphericalDirection(sintheta, costheta, phi);
    if (!Vector.SameHemisphere(wo, wh)) {
      wh = -wh;
    }

    // Compute incident direction by reflecting about [wh]
    wi.copy(-wo + wh * 2.0 * Vector.Dot(wo, wh));

    // Compute PDF for [wi] from Blinn distribution
    double blinn_pdf = ((exponent + 1.0) * Math.pow(costheta, exponent)) /
                       (2.0 * Math.PI * 4.0 * Vector.Dot(wo, wh));

    if (Vector.Dot(wo, wh) <= 0.0) {
      blinn_pdf = 0.0;
    }

    return blinn_pdf;
  }

  double pdf(Vector wo, Vector wi) {
    Vector wh = Vector.Normalize(wo + wi);
    double costheta = Vector.AbsCosTheta(wh);
    // Compute PDF for [wi] from Blinn distribution
    double blinn_pdf = ((exponent + 1.0) * Math.pow(costheta, exponent)) /
                      (2.0 * Math.PI * 4.0 * Vector.Dot(wo, wh));
    if (Vector.Dot(wo, wh) <= 0.0) {
      blinn_pdf = 0.0;
    }
    return blinn_pdf;
  }

  double exponent;
}
