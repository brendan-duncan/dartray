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

/**
 * Defines a reflection model for dielectric media
 * (objects that don't conduct electricity, like glass).
 *
 * The parameters {eta_i} and {eta_t} define the refraction on the
 * outside and the inside of the surface respectively.
 */
class FresnelDielectric extends Fresnel {
  FresnelDielectric(this.eta_i, this.eta_t);

  Spectrum evaluate(double cosi) {
    // Compute Fresnel reflectance for dielectric
    cosi = cosi.clamp(-1.0, 1.0);

    // Compute indices of refraction for dielectric
    bool entering = cosi > 0.0;
    double ei = eta_i;
    double et = eta_t;
    if (!entering) {
      double t = ei;
      ei = et;
      et = t;
    }

    // Compute _sint_ using Snell's law
    double sint = ei / et * Math.sqrt(Math.max(0.0, 1.0 - cosi * cosi));
    if (sint >= 1.0) {
      // Handle total internal reflection
      return new Spectrum(1.0);
    }

    double cost = Math.sqrt(Math.max(0.0, 1.0 - sint * sint));
    cosi = cosi.abs();

    double Rparl = ((et * cosi) - (ei * cost)) /
                   ((et * cosi) + (ei * cost));

    double Rperp = ((ei * cosi) - (et * cost)) /
                   ((ei * cosi) + (et * cost));

    return new Spectrum((Rparl * Rparl + Rperp * Rperp) / 2.0);
  }

  double eta_i;
  double eta_t;
}
