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

class SphericalMapping2D extends TextureMapping2D {
  SphericalMapping2D(this.worldToTexture);

  void map(DifferentialGeometry dg, List<double> s, List<double> t,
             List<double> dsdx, List<double> dtdx, List<double> dsdy,
             List<double> dtdy) {
    _sphere(dg.p, s, t);
    // Compute texture coordinate differentials for sphere $(u,v)$ mapping
    List<double> sx = [0.0];
    List<double> tx = [0.0];
    List<double> sy = [0.0];
    List<double> ty = [0.0];

    const double delta = 0.1;

    _sphere(dg.p + dg.dpdx * delta, sx, tx);
    dsdx[0] = (sx[0] - s[0]) / delta;
    dtdx[0] = (tx[0] - t[0]) / delta;

    if (dtdx[0] > 0.5) {
      dtdx[0] = 1.0 - dtdx[0];
    } else if (dtdx[0] < -0.5) {
      dtdx[0] = -(dtdx[0] + 1.0);
    }

    _sphere(dg.p + dg.dpdy * delta, sy, ty);
    dsdy[0] = (sy[0] - s[0]) / delta;
    dtdy[0] = (ty[0] - t[0]) / delta;

    if (dtdy[0] > 0.5) {
      dtdy[0] = 1.0 - dtdy[0];
    } else if (dtdy[0] < -0.5) {
      dtdy[0] = -(dtdy[0] + 1.0);
    }
  }

  void _sphere(Point p, List<double> s, List<double> t) {
    Vector vec = Vector.Normalize(worldToTexture.transformPoint(p));
    double theta = Vector.SphericalTheta(vec);
    double phi = Vector.SphericalPhi(vec);
    s[0] = theta * INV_PI;
    t[0] = phi * INV_TWOPI;
  }

  Transform worldToTexture;
}
