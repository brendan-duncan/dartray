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

abstract class VolumeRegion {
  BBox worldBound();

  bool intersectP(Ray ray, List<double> t0, List<double> t1);

  Spectrum sigma_a(Point p, Vector w, double time);

  Spectrum sigma_s(Point p, Vector w, double time);

  Spectrum Lve(Point p, Vector w, double time);

  double p(Point p, Vector w, Vector wp, double time);

  Spectrum sigma_t(Point p, Vector wo, double time) {
    return sigma_a(p, wo, time) + sigma_s(p, wo, time);
  }

  Spectrum tau(Ray ray, [double step = 1.0, double offset = 0.5]);
}
