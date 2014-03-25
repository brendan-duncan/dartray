/****************************************************************************
 *  Copyright (C) 2014 by Brendan Duncan.                                   *
 *                                                                          *
 *  This file is part of DartRay.                                           *
 *                                                                          *
 *  Licensed under the Apache License, Version 2.0 (the "License");         *
 *  you may not use this file except in compliance with the License.        *
 *  You may obtain a copy of the License at                                 *
 *                                                                          *
 *  http://www.apache.org/licenses/LICENSE-2.0                              *
 *                                                                          *
 *  Unless required by applicable law or agreed to in writing, software     *
 *  distributed under the License is distributed on an "AS IS" BASIS,       *
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 *  See the License for the specific language governing permissions and     *
 *  limitations under the License.                                          *
 *                                                                          *
 *   This project is based on PBRT v2 ; see http://www.pbrt.org             *
 *   pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.*
 ****************************************************************************/
part of core;

class SpecularReflection extends BxDF {
  SpecularReflection(RGBColor r, Fresnel f) :
    super(BSDF_REFLECTION | BSDF_SPECULAR),
    R = new RGBColor.from(r),
    fresnel = f;

  RGBColor f(Vector wo, Vector wi) {
    return new RGBColor(0.0);
  }

  RGBColor sample_f(Vector wo, Vector wi, double u1, double u2,
                       List<double> pdf) {
    // Compute perfect specular reflection direction
    wi.x = -wo.x;
    wi.y = -wo.y;
    wi.z = wo.z;
    pdf[0] = 1.0;
    return (fresnel.evaluate(Vector.CosTheta(wo)) * R).invScale(Vector.AbsCosTheta(wi));
  }

  double pdf(Vector wo, Vector wi) {
    return 0.0;
  }

  RGBColor R;
  Fresnel fresnel;
}
