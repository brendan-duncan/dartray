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

class CheckerboardTexture extends Texture {
  static const int AA_NONE = 0;
  static const int AA_CLOSED_FORM = 1;

  CheckerboardTexture(this.mapping, this.tex1, this.tex2, this.aaMethod);

  static Texture _Create(Transform tex2world, TextureParams tp,
                         Texture tex1, Texture tex2) {
    int dim = tp.findInt("dimension", 2);
    if (dim != 2 && dim != 3) {
      LogWarning("$dim dimensional checkerboard texture not supported");
      return null;
    }

    if (dim == 2) {
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
                                  tp.findFloat("udelta", 0.0),
                                  tp.findFloat("vdelta", 0.0));
      } else {
        LogWarning("2D texture mapping \"$type\" unknown");
        map = new UVMapping2D();
      }

      String aa = tp.findString("aamode", "closedform");
      int aaMethod = AA_CLOSED_FORM;
      if (aa == "none") {
        aaMethod = AA_NONE;
      } else if (aa == "closedform") {
        aaMethod = AA_CLOSED_FORM;
      } else {
        LogWarning("Antialiasing mode \"$aa\" not understood by "
                 "Checkerboard2DTexture; using \"closedform\"");
        aaMethod = AA_CLOSED_FORM;
      }
      return new CheckerboardTexture(map, tex1, tex2, aaMethod);
    } else {
      // Initialize 3D texture mapping _map_ from _tp_
      TextureMapping3D map = new IdentityMapping3D(tex2world);
      return new Checkerboard3DTexture(map, tex1, tex2);
    }
  }

  static Texture CreateFloat(Transform tex2world, TextureParams tp) {
    Texture tex1 = tp.getFloatTexture("tex1", 1.0);
    Texture tex2 = tp.getFloatTexture("tex2", 0.0);
    return _Create(tex2world, tp, tex1, tex2);
  }

  static Texture CreateSpectrum(Transform tex2world, TextureParams tp) {
    Texture tex1 = tp.getSpectrumTexture("tex1", new Spectrum(1.0));
    Texture tex2 = tp.getSpectrumTexture("tex2", new Spectrum(0.0));
    return _Create(tex2world, tp, tex1, tex2);
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
    double dsdx = _dsdx[0];
    double dsdy = _dsdy[0];
    double dtdx = _dtdx[0];
    double dtdy = _dtdy[0];

    if (aaMethod == AA_NONE) {
      // Point sample _Checkerboard2DTexture_
      if ((s.floor() + t.floor()) % 2 == 0) {
        return tex1.evaluate(dg);
      }
      return tex2.evaluate(dg);
    } else {
      // Compute closed-form box-filtered _Checkerboard2DTexture_ value

      // Evaluate single check if filter is entirely inside one of them
      double ds = max(dsdx.abs(), dsdy.abs());
      double dt = max(dtdx.abs(), dtdy.abs());
      double s0 = s - ds;
      double s1 = s + ds;
      double t0 = t - dt;
      double t1 = t + dt;
      if (s0.floor() == s1.floor() && t0.floor() == t1.floor()) {
        // Point sample _Checkerboard2DTexture_
        if ((s.floor() + t.floor()) % 2 == 0) {
          return tex1.evaluate(dg);
        }
        return tex2.evaluate(dg);
      }

      // Apply box filter to checkerboard region
      BUMPINT(x) => ((x / 2).floor() +
                     2.0 * max((x / 2) - (x / 2).floor() - 0.5, 0.0));
      double sint = (BUMPINT(s1) - BUMPINT(s0)) / (2.0 * ds);
      double tint = (BUMPINT(t1) - BUMPINT(t0)) / (2.0 * dt);
      double area2 = sint + tint - 2.0 * sint * tint;
      if (ds > 1.0 || dt > 1.0) {
        area2 = 0.5;
      }

      return tex1.evaluate(dg) * (1.0 - area2) + tex2.evaluate(dg) * area2;
    }
  }

  TextureMapping2D mapping;
  Texture tex1;
  Texture tex2;
  int aaMethod;
}
