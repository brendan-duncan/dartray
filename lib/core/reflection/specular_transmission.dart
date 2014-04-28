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

class SpecularTransmission extends BxDF {
  SpecularTransmission(Spectrum t, double ei, double et)
      : fresnel = new FresnelDielectric(ei, et),
        super(BSDF_TRANSMISSION | BSDF_SPECULAR) {
    T = t;
    etai = ei;
    etat = et;
  }

  Spectrum f(Vector wo, Vector wi) {
    return new Spectrum(0.0);
  }

  Spectrum sample_f(Vector wo, Vector wi, double u1, double u2,
                    List<double> pdf) {
    // Figure out which [eta] is incident and which is transmitted
    bool entering = Vector.CosTheta(wo) > 0.0;
    double ei = etai;
    double et = etat;
    if (!entering) {
      double t = ei;
      ei = et;
      et = t;
    }

    // Compute transmitted ray direction
    double sini2 = Vector.SinTheta2(wo);
    double eta = ei / et;
    double sint2 = eta * eta * sini2;

    // Handle total internal reflection for transmission
    if (sint2 >= 1.0) {
      return new Spectrum(0.0);
    }

    double cost = Math.sqrt(Math.max(0.0, 1.0 - sint2));
    if (entering) {
      cost = -cost;
    }
    double sintOverSini = eta;
    wi.x = sintOverSini * -wo.x;
    wi.y = sintOverSini * -wo.y;
    wi.z = cost;
    pdf[0] = 1.0;
    Spectrum F = fresnel.evaluate(Vector.CosTheta(wo));

    return ((Spectrum.ONE - F) * T) / Vector.AbsCosTheta(wi);
  }

  double pdf(Vector wo, Vector wi) {
    return 0.0;
  }

  Spectrum T;
  double etai, etat;
  FresnelDielectric fresnel;
}
