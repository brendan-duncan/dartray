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

class FBmTexture extends Texture {
  FBmTexture(this.octaves, this.omega, this.mapping, this.spectrum);

  static FBmTexture CreateFloat(Transform tex2world, TextureParams tp) {
    // Initialize 3D texture mapping _map_ from _tp_
    TextureMapping3D map = new IdentityMapping3D(tex2world);
    return new FBmTexture(tp.findInt('octaves', 8),
                          tp.findFloat('roughness', 0.5),
                          map, false);
  }

  static FBmTexture CreateSpectrum(Transform tex2world, TextureParams tp) {
    // Initialize 3D texture mapping _map_ from _tp_
    TextureMapping3D map = new IdentityMapping3D(tex2world);
    return new FBmTexture(tp.findInt('octaves', 8),
                          tp.findFloat('roughness', 0.5),
                          map, true);
  }

  evaluate(DifferentialGeometry dg) {
    Vector dpdx = new Vector();
    Vector dpdy = new Vector();
    Point P = mapping.map(dg, dpdx, dpdy);
    double n = FBm(P, dpdx, dpdy, omega, octaves);
    return spectrum ? new Spectrum(n) : n;
  }

  double omega;
  int octaves;
  bool spectrum;
  TextureMapping3D mapping;
}
