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
part of textures;

class MarbleTexture extends Texture {
  MarbleTexture(this.octaves, this.omega, this.scale, this.variation,
                this.mapping);

  static MarbleTexture CreateFloat(Transform tex2world, TextureParams tp) {
    return null;
  }

  static MarbleTexture CreateSpectrum(Transform tex2world, TextureParams tp) {
    // Initialize 3D texture mapping _map_ from _tp_
    TextureMapping3D map = new IdentityMapping3D(tex2world);
    return new MarbleTexture(tp.findInt('octaves', 8),
                          tp.findFloat('roughness', 0.5),
                          tp.findFloat('scale', 1.0),
                          tp.findFloat('variation', 0.2),
                          map);
  }

  evaluate(DifferentialGeometry dg) {
    Vector dpdx = new Vector();
    Vector dpdy = new Vector();
    Point P = mapping.map(dg, dpdx, dpdy);
    P *= scale;
    double marble = P.y + variation *
                   FBm(P, dpdx * scale, dpdy * scale, omega, octaves);
    double t = 0.5 + 0.5 * sin(marble);

    // Evaluate marble spline at _t_
    const List<double> c = const [ 0.58, 0.58, 0.6,
                                   0.58, 0.58, 0.6,
                                   0.58, 0.58, 0.6,
                                   0.5, 0.5, 0.5,
                                   0.6, 0.59, 0.58,
                                   0.58, 0.58, 0.6,
                                   0.58, 0.58, 0.6,
                                   0.2, 0.2, 0.33,
                                   0.58, 0.58, 0.6];
    const int NC = 9;
    const int NSEG = NC - 3;
    int first = (t * NSEG).floor();
    t = (t * NSEG - first);
    int ci = first * 3;
    Spectrum c0 = new Spectrum.rgb(c[ci], c[ci + 1], c[ci + 2]);
    Spectrum c1 = new Spectrum.rgb(c[ci + 3], c[ci + 4], c[ci + 5]);
    Spectrum c2 = new Spectrum.rgb(c[ci + 6], c[ci + 7], c[ci + 8]);
    Spectrum c3 = new Spectrum.rgb(c[ci + 9], c[ci + 10], c[ci + 11]);

    // Bezier spline evaluated with de Castilejau's algorithm
    Spectrum s0 = c0 * (1.0 - t) + c1 * t;
    Spectrum s1 = c1 * (1.0 - t) + c2 * t;
    Spectrum s2 = c2 * (1.0 - t) + c3 * t;
    s0 = s0 * (1.0 - t) + s1 * t;
    s1 = s1 * (1.0 - t) + s2 * t;
    // Extra scale of 1.5 to increase variation among colors
    return (s0 * (1.0 - t) + s1 * t) * 1.5;
  }

  int octaves;
  double omega, scale, variation;
  TextureMapping3D mapping;
}
