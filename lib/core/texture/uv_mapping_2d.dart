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

class UVMapping2D extends TextureMapping2D {
  UVMapping2D([this.su = 1.0, this.sv = 1.0, this.du = 0.0, this.dv = 0.0]);

  void map(DifferentialGeometry dg, List<double> s, List<double> t,
           List<double> dsdx, List<double> dtdx, List<double> dsdy,
           List<double> dtdy) {
    s[0] = su * dg.u + du;
    t[0] = sv * dg.v + dv;
    // Compute texture differentials for 2D identity mapping
    dsdx[0] = su * dg.dudx;
    dtdx[0] = sv * dg.dvdx;
    dsdy[0] = su * dg.dudy;
    dtdy[0] = sv * dg.dvdy;
  }

  double su;
  double sv;
  double du;
  double dv;
}
