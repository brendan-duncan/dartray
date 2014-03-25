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

class FresnelConductor extends Fresnel {
  FresnelConductor(RGBColor e, RGBColor kk) :
    eta = new RGBColor.from(e),
    k = new RGBColor.from(kk);

  RGBColor evaluate(double cosi) {
    RGBColor one = new RGBColor(1.0);
    RGBColor cosSqr = new RGBColor(cosi * cosi);

    RGBColor tmp = (eta * eta + k * k) * (cosi * cosi);
    RGBColor Rparl2 = (tmp - (eta * (2.0 * cosi)) + one) /
                         (tmp + (eta * (2.0 * cosi)) + one);
    RGBColor tmp_f = eta * eta + k * k;
    RGBColor Rperp2 = (tmp_f - (eta * (2.0 * cosi)) + cosSqr) /
                         (tmp_f + (eta * (2.0 * cosi)) + cosSqr);

    return (Rparl2 + Rperp2) / 2.0;
  }

  RGBColor eta;
  RGBColor k;
}
