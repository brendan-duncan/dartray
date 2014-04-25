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

class ImageTexture extends Texture {
  ImageTexture(this.mapping, String filename, bool doTri,
               double maxAniso, int wrap, double scale, double gamma,
               bool spectrum) {
    if (filename.isNotEmpty) {
      Completer completer = new Completer();
      ResourceManager.RequestImage(filename, completer.future)
        .then((SpectrumImage img) {
          String name = MIPMap.GetTextureName(filename, doTri: doTri,
                                              maxAniso: maxAniso,
                                              wrap: wrap, scale: scale,
                                              gamma: gamma,
                                              spectrum: spectrum);
          LogDebug('TEXTURE $name');

          if (ResourceManager.HasTexture(name)) {
            mipmap = ResourceManager.GetTexture(name);
            completer.complete();
            return;
          }

          if (img != null) {
            if (!spectrum) {
              img = img.convert(SpectrumImage.FLOAT);
            } else {
              img = new SpectrumImage.from(img);
            }

            if (scale != 1.0 || gamma != 1.0) {
              for (int i = 0, len = img.data.length; i < len; ++i) {
                img.data[i] = pow(img.data[i] * scale, gamma);
              }
            }

            mipmap = new MIPMap.texture(img, filename, doTri, maxAniso, wrap);
            ResourceManager.AddTexture(name, mipmap);
          }

          // Let the renderer know we're done processing the resource.
          completer.complete();
        });
    }

    double v = pow(scale, gamma);
    SpectrumImage img = new SpectrumImage(1, 1, spectrum ?
                                          SpectrumImage.SPECTRUM :
                                          SpectrumImage.FLOAT);
    if (spectrum) {
      img[0] = new Spectrum(v);
    } else {
      img[0] = v;
    }

    mipmap = new MIPMap.texture(img, '');
  }

  evaluate(DifferentialGeometry dg) {
    List<double> s = [0.0];
    List<double> t = [0.0];
    List<double> dsdx = [0.0];
    List<double> dtdx = [0.0];
    List<double> dsdy = [0.0];
    List<double> dtdy = [0.0];
    mapping.map(dg, s, t, dsdx, dtdx, dsdy, dtdy);
    var v = mipmap.lookup2(s[0], t[0], dsdx[0], dtdx[0], dsdy[0], dtdy[0]);
    return v;
  }

  MIPMap mipmap;
  TextureMapping2D mapping;


  static ImageTexture CreateFloat(Transform tex2world, TextureParams tp) {
    // Initialize 2D texture mapping _map_ from _tp_
    TextureMapping2D map;
    String type = tp.findString('mapping', 'uv');
    if (type == 'uv') {
      double su = tp.findFloat('uscale', 1.0);
      double sv = tp.findFloat('vscale', 1.0);
      double du = tp.findFloat('udelta', 0.0);
      double dv = tp.findFloat('vdelta', 0.0);
      map = new UVMapping2D(su, sv, du, dv);
    } else if (type == 'spherical') {
      map = new SphericalMapping2D(Transform.Inverse(tex2world));
    } else if (type == 'cylindrical') {
      map = new CylindricalMapping2D(Transform.Inverse(tex2world));
    } else if (type == 'planar') {
      map = new PlanarMapping2D(tp.findVector('v1', new Vector(1.0, 0.0, 0.0)),
                                tp.findVector('v2', new Vector(0.0, 1.0, 0.0)),
                                tp.findFloat('udelta', 0.0),
                                tp.findFloat('vdelta', 0.0));
    } else {
      LogError('2D texture mapping \'$type\' unknown');
      map = new UVMapping2D();
    }

    // Initialize _ImageTexture_ parameters
    double maxAniso = tp.findFloat('maxanisotropy', 8.0);
    bool trilerp = tp.findBool('trilinear', false);
    String wrap = tp.findString('wrap', 'repeat');
    int wrapMode = MIPMap.TEXTURE_REPEAT;
    if (wrap == 'black') {
      wrapMode = MIPMap.TEXTURE_BLACK;
    } else if (wrap == 'clamp') {
      wrapMode = MIPMap.TEXTURE_CLAMP;
    }
    double scale = tp.findFloat('scale', 1.0);
    double gamma = tp.findFloat('gamma', 1.0);
    return new ImageTexture(map, tp.findFilename('filename'),
          trilerp, maxAniso, wrapMode, scale, gamma, false);
  }

  static ImageTexture CreateSpectrum(Transform tex2world, TextureParams tp) {
    // Initialize 2D texture mapping _map_ from _tp_
    TextureMapping2D map;
    String type = tp.findString('mapping', 'uv');
    if (type == 'uv') {
      double su = tp.findFloat('uscale', 1.0);
      double sv = tp.findFloat('vscale', 1.0);
      double du = tp.findFloat('udelta', 0.0);
      double dv = tp.findFloat('vdelta', 0.0);
      map = new UVMapping2D(su, sv, du, dv);
    } else if (type == 'spherical') {
      map = new SphericalMapping2D(Transform.Inverse(tex2world));
    } else if (type == 'cylindrical') {
      map = new CylindricalMapping2D(Transform.Inverse(tex2world));
    } else if (type == 'planar') {
      map = new PlanarMapping2D(tp.findVector('v1', new Vector(1.0, 0.0, 0.0)),
                                tp.findVector('v2', new Vector(0.0, 1.0, 0.0)),
                                tp.findFloat('udelta', 0.0),
                                tp.findFloat('vdelta', 0.0));
    } else {
      LogError('2D texture mapping \'$type\' unknown');
      map = new UVMapping2D();
    }

    // Initialize _ImageTexture_ parameters
    double maxAniso = tp.findFloat('maxanisotropy', 8.0);
    bool trilerp = tp.findBool('trilinear', false);
    String wrap = tp.findString('wrap', 'repeat');
    int wrapMode = MIPMap.TEXTURE_REPEAT;
    if (wrap == 'black') {
      wrapMode = MIPMap.TEXTURE_BLACK;
    } else if (wrap == 'clamp') {
      wrapMode = MIPMap.TEXTURE_CLAMP;
    }

    double scale = tp.findFloat('scale', 1.0);
    double gamma = tp.findFloat('gamma', 1.0);

    return new ImageTexture(map, tp.findFilename('filename'),
                            trilerp, maxAniso, wrapMode, scale, gamma, true);
  }
}
