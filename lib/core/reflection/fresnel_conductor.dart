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

class FresnelConductor extends Fresnel {
  FresnelConductor(Spectrum e, Spectrum kk)
      : eta = new Spectrum.from(e),
        k = new Spectrum.from(kk);

  Spectrum evaluate(double cosi) {
    cosi = cosi.abs();

    Spectrum cosSqr = new Spectrum(cosi * cosi);

    Spectrum tmp = (eta * eta + k * k) * (cosi * cosi);
    Spectrum r1 = (tmp - (eta * (2.0 * cosi)) + ONE);
    Spectrum r2 = (tmp + (eta * (2.0 * cosi)) + ONE);
    Spectrum Rparl2 = r1 / r2;

    Spectrum tmp_f = eta * eta + k * k;
    r1 = (tmp_f - (eta * (2.0 * cosi)) + cosSqr);
    r2 = (tmp_f + (eta * (2.0 * cosi)) + cosSqr);
    Spectrum Rperp2 = r1 / r2;

    return (Rparl2 + Rperp2) / 2.0;
  }

  static final Spectrum ONE = new Spectrum(1.0);
  Spectrum eta;
  Spectrum k;
}
