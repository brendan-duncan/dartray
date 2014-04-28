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

class PlanarMapping2D extends TextureMapping2D {
  PlanarMapping2D(Vector vv1, Vector vv2, [this.ds = 0.0, this.dt = 0.0])
      : vs = new Vector.from(vv1),
        vt = new Vector.from(vv2);

  void map(DifferentialGeometry dg, List<double> s, List<double> t,
           List<double> dsdx, List<double> dtdx, List<double> dsdy,
           List<double> dtdy) {
    Vector vec = dg.p;
    s[0] = ds + Vector.Dot(vec, vs);
    t[0] = dt + Vector.Dot(vec, vt);
    dsdx[0] = Vector.Dot(dg.dpdx, vs);
    dtdx[0] = Vector.Dot(dg.dpdx, vt);
    dsdy[0] = Vector.Dot(dg.dpdy, vs);
    dtdy[0] = Vector.Dot(dg.dpdy, vt);
  }

  Vector vs;
  Vector vt;
  double ds;
  double dt;
}
