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

class UVTexture extends Texture {
  UVTexture(this.mapping);

  static UVTexture CreateFloat(Transform tex2world, TextureParams tp) {
    return null;
  }

  static UVTexture CreateSpectrum(Transform tex2world, TextureParams tp) {
    // Initialize 2D texture mapping _map_ from _tp_
    TextureMapping2D map;
    String type = tp.findString("mapping", "uv");
    if (type == "uv") {
      double su = tp.findFloat("uscale", 1.0);
      double sv = tp.findFloat("vscale", 1.0);
      double du = tp.findFloat("udelta", 0.0);
      double dv = tp.findFloat("vdelta", 0.0);
      map = new UVMapping2D(su, sv, du, dv);
    } else if (type == "spherical") {
      map = new SphericalMapping2D(Transform.Inverse(tex2world));
    } else if (type == "cylindrical") {
      map = new CylindricalMapping2D(Transform.Inverse(tex2world));
    } else if (type == "planar") {
      map = new PlanarMapping2D(tp.findVector("v1", new Vector(1.0, 0.0, 0.0)),
          tp.findVector("v2", new Vector(0.0, 1.0, 0.0)),
          tp.findFloat("udelta", 0.0), tp.findFloat("vdelta", 0.0));
    } else {
      LogError("2D texture mapping \"$type\" unknown");
      map = new UVMapping2D();
    }

    return new UVTexture(map);
  }

  evaluate(DifferentialGeometry dg) {
    List<double> _s = [0.0],
                 _t = [0.0],
                 _dsdx = [0.0],
                 _dtdx = [0.0],
                 _dsdy = [0.0],
                 _dtdy = [0.0];
    mapping.map(dg, _s, _t, _dsdx, _dtdx, _dsdy, _dtdy);
    double s = _s[0];
    double t = _t[0];
    return new Spectrum.rgb(s - s.floor(), t - t.floor(), 0.0);
  }

  TextureMapping2D mapping;
}
