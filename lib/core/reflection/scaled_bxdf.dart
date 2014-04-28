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

class ScaledBxDF extends BxDF {
  ScaledBxDF(BxDF b, Spectrum sc)
      : bxdf = b,
        s = new Spectrum.from(sc),
        super(b.type);

  Spectrum rho(Vector w, int nSamples, List<double> samples) {
    return s * bxdf.rho(w, nSamples, samples);
  }

  Spectrum rho2(int nSamples, List<double> samples1, List<double> samples2) {
    return s * bxdf.rho2(nSamples, samples1, samples2);
  }

  Spectrum f(Vector wo, Vector wi) {
    return s * bxdf.f(wo, wi);
  }

  Spectrum sample_f(Vector wo, Vector wi, double u1, double u2,
                    List<double> pdf) {
    Spectrum f = bxdf.sample_f(wo, wi, u1, u2, pdf);
    return s * f;
  }

  BxDF bxdf;
  Spectrum s;
}
