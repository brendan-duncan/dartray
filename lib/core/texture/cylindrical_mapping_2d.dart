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

class CylindricalMapping2D extends TextureMapping2D {
  CylindricalMapping2D(this.worldToTexture);

  void map(DifferentialGeometry dg, List<double> s, List<double> t,
           List<double> dsdx, List<double> dtdx, List<double> dsdy,
           List<double> dtdy) {
    _cylinder(dg.p, s, t);
    // Compute texture coordinate differentials for cylinder (u,v) mapping
    List<double> sx = [0.0],
                 tx = [0.0],
                 sy = [0.0],
                 ty = [0.0];

    double delta = 0.01;
    _cylinder(dg.p + dg.dpdx * delta, sx, tx);
    dsdx[0] = (sx[0] - s[0]) / delta;
    dtdx[0] = (tx[0] - t[0]) / delta;
    if (dtdx[0] > 0.5) {
      dtdx[0] = 1.0 - dtdx[0];
    } else if (dtdx[0] < -0.5) {
      dtdx[0] = -(dtdx[0] + 1.0);
    }

    _cylinder(dg.p + dg.dpdy * delta, sy, ty);
    dsdy[0] = (sy[0] - s[0]) / delta;
    dtdy[0] = (ty[0] - t[0]) / delta;
    if (dtdy[0] > 0.5) {
      dtdy[0] = 1.0 - dtdy[0];
    } else if (dtdy[0] < -0.5) {
      dtdy[0] = -(dtdy[0] + 1.0);
    }
  }

  void _cylinder(Point p, List<double> s, List<double> t) {
    Vector vec = Vector.Normalize(worldToTexture.transformPoint(p));
    s[0] = (Math.PI + Math.atan2(vec.y, vec.x)) / (2.0 * Math.PI);
    t[0] = vec.z;
  }

  Transform worldToTexture;
}
