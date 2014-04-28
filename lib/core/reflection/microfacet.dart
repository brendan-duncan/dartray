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

class Microfacet extends BxDF {
  Microfacet(this.R, this.fresnel, this.distribution)
      : super(BSDF_REFLECTION | BSDF_GLOSSY);

  Spectrum f(Vector wo, Vector wi) {
    double cosThetaO = Vector.AbsCosTheta(wo);
    double cosThetaI = Vector.AbsCosTheta(wi);
    if (cosThetaI == 0.0 || cosThetaO == 0.0) {
      return new Spectrum(0.0);
    }

    Vector wh = wi + wo;
    if (wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0) {
      return new Spectrum(0.0);
    }

    wh = Vector.Normalize(wh);

    double cosThetaH = Vector.Dot(wi, wh);
    Spectrum F = fresnel.evaluate(cosThetaH);

    return R * (distribution.d(wh) * g(wo, wi, wh)) *
           F / (4.0 * cosThetaI * cosThetaO);
  }

  double g(Vector wo, Vector wi, Vector wh) {
    double NdotWh = Vector.AbsCosTheta(wh);
    double NdotWo = Vector.AbsCosTheta(wo);
    double NdotWi = Vector.AbsCosTheta(wi);
    double WOdotWh = Vector.AbsDot(wo, wh);
    return Math.min(1.0, Math.min((2.0 * NdotWh * NdotWo / WOdotWh),
                        (2.0 * NdotWh * NdotWi / WOdotWh)));
  }

  Spectrum sample_f(Vector wo, Vector wi, double u1, double u2,
                    List<double> pdf) {
    pdf[0] = distribution.sample_f(wo, wi, u1, u2);
    if (!Vector.SameHemisphere(wo, wi)) {
      return new Spectrum(0.0);
    }
    return f(wo, wi);
  }

  double pdf(Vector wo, Vector wi) {
    if (!Vector.SameHemisphere(wo, wi)) {
      return 0.0;
    }
    return distribution.pdf(wo, wi);
  }

  Spectrum R;
  MicrofacetDistribution distribution;
  Fresnel fresnel;
}
