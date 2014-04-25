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


class TextureParams {
  TextureParams(this.geomParams, this.materialParams, this.floatTextures,
                this.spectrumTextures);

  Texture getSpectrumTexture(String n, Spectrum def) {
    String name = geomParams.findTexture(n);
    if (name.isEmpty) {
      name = materialParams.findTexture(n);
    }
    if (name.isNotEmpty) {
      if (spectrumTextures.containsKey(name)) {
          return spectrumTextures[name];
      } else {
          LogWarning('Couldn\'t find spectrum texture named \'$name\' '
                     'for parameter \'$n\'');
      }
    }

    Spectrum val = geomParams.findOneSpectrum(n,
                                materialParams.findOneSpectrum(n, def));

    return new ConstantTexture(val);
  }

  Texture getFloatTexture(String n, double def) {
    String name = geomParams.findTexture(n);
    if (name == '') {
      name = materialParams.findTexture(n);
    }
    if (name != '') {
      if (floatTextures.containsKey(name)) {
        return floatTextures[name];
      } else {
        LogWarning('Couldn\'t find float texture named \'$name\' '
                   'for parameter \'$n\'');
      }
    }

    double val = geomParams.findOneFloat(n,
        materialParams.findOneFloat(n, def));

    return new ConstantTexture(val);
  }

  Texture getFloatTextureOrNull(String n) {
    String name = geomParams.findTexture(n);
    if (name == '') {
      name = materialParams.findTexture(n);
    }
    if (name == '') {
      return null;
    }

    if (floatTextures.containsKey(name)) {
      return floatTextures[name];
    } else {
      LogWarning('Couldn\'t find float texture named \'$name\' '
                 'for parameter \'$n\'');
      return null;
    }
  }

  double findFloat(String n, double d) {
    return geomParams.findOneFloat(n, materialParams.findOneFloat(n, d));
  }

  String findString(String n, [String d = '']) {
    return geomParams.findOneString(n, materialParams.findOneString(n, d));
  }

  String findFilename(String n, [String d = '']) {
    return geomParams.findOneFilename(n, materialParams.findOneFilename(n, d));
  }

  int findInt(String n, int d) {
    return geomParams.findOneInt(n, materialParams.findOneInt(n, d));
  }

  bool findBool(String n, bool d) {
    return geomParams.findOneBool(n, materialParams.findOneBool(n, d));
  }

  Point findPoint(String n, Point d) {
    return geomParams.findOnePoint(n, materialParams.findOnePoint(n, d));
  }

  Vector findVector(String n, Vector d) {
    return geomParams.findOneVector(n, materialParams.findOneVector(n, d));
  }

  Normal findNormal(String n, Normal d) {
    return geomParams.findOneNormal(n, materialParams.findOneNormal(n, d));
  }

  Spectrum findSpectrum(String n, Spectrum d) {
    return geomParams.findOneSpectrum(n, materialParams.findOneSpectrum(n, d));
  }

  void reportUnused() {
    geomParams.reportUnused();
    materialParams.reportUnused();
  }

  ParamSet getGeomParams() {
    return geomParams;
  }

  ParamSet getMaterialParams() {
    return materialParams;
  }

  Map<String, Texture> floatTextures;
  Map<String, Texture> spectrumTextures;

  ParamSet geomParams;
  ParamSet materialParams;
}
