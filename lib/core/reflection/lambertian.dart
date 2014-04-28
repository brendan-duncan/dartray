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
 * The [Lambertian] BRDF models a perfectly diffuse surface that scatters
 * incident illumination equally in all directions.
 *
 * Altough this reflection model isn't phsyically plausible, it is a good
 * approximation to many real-world surfaces such as matte paint.
 */
class Lambertian extends BxDF {
  Lambertian(Spectrum reflectance)
      : R = new Spectrum.from(reflectance),
        super(BSDF_REFLECTION | BSDF_DIFFUSE);

  Spectrum f(Vector wo, Vector wi) {
    return R * INV_PI;
  }

  Spectrum rho(Vector wo, int nSamples, List<double> samples) {
    return R;
  }

  Spectrum rho2(int nSamples, List<double> samples1, List<double> samples2) {
    return R;
  }

  Spectrum R;
}
