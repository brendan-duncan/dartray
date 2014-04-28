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

class OrenNayar extends BxDF {
  OrenNayar(Spectrum reflectance, double sig)
      : R = new Spectrum.from(reflectance),
        super(BSDF_REFLECTION | BSDF_DIFFUSE) {
    double sigma = Radians(sig);
    double sigma2 = sigma * sigma;
    A = 1.0 - (sigma2 / (2.0 * (sigma2 + 0.33)));
    B = 0.45 * sigma2 / (sigma2 + 0.09);
  }

  Spectrum f(Vector wo, Vector wi) {
    double sinthetai = Vector.SinTheta(wi);
    double sinthetao = Vector.SinTheta(wo);
    // Compute cosine term of Oren-Nayar model
    double maxcos = 0.0;
    if (sinthetai > 1e-4 && sinthetao > 1e-4) {
      double sinphii = Vector.SinPhi(wi);
      double cosphii = Vector.CosPhi(wi);
      double sinphio = Vector.SinPhi(wo);
      double cosphio = Vector.CosPhi(wo);
      double dcos = cosphii * cosphio + sinphii * sinphio;
      maxcos = Math.max(0.0, dcos);
    }

    // Compute sine and tangent terms of Oren-Nayar model
    double sinalpha, tanbeta;
    if (Vector.AbsCosTheta(wi) > Vector.AbsCosTheta(wo)) {
      sinalpha = sinthetao;
      tanbeta = sinthetai / Vector.AbsCosTheta(wi);
    } else {
      sinalpha = sinthetai;
      tanbeta = sinthetao / Vector.AbsCosTheta(wo);
    }

    return R * (INV_PI * (A + B * maxcos * sinalpha * tanbeta));
  }

  Spectrum R;
  double A;
  double B;
}
