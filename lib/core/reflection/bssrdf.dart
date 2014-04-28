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
 * Bidirectional surface scattering reflectance distribution function describes
 * the relation between outgoing radiance and the incident flux, including the
 * phenomena like subsurface scattering (SSS). The BSSRDF describes how light
 * is transported between any two rays that hit a surface.
 */
class BSSRDF {
  final double eta;
  final Spectrum sigma_a;
  final Spectrum sigma_prime_s;

  BSSRDF(Spectrum sigma_a, Spectrum sigma_prime_s, this.eta)
      : this.sigma_a = new Spectrum.from(sigma_a),
        this.sigma_prime_s = new Spectrum.from(sigma_prime_s);
}
