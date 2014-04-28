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

class ParamSet {
  ParamSet();

  ParamSet.from(ParamSet other)
      : bools = new List.from(other.bools),
        ints = new List.from(other.ints),
        floats = new List.from(other.floats),
        points = new List.from(other.points),
        vectors = new List.from(other.vectors),
        normals = new List.from(other.normals),
        spectra = new List.from(other.spectra),
        strings = new List.from(other.strings),
        textures = new List.from(other.textures);

  ParamSet.fromJson(Map json) {
    for (String type_name in json.keys) {
      var value = json[type_name];
      List tk = type_name.split(' ');
      if (tk.length != 2) {
        LogError('Invalid parameter declaration: \'$type_name\'. '
                 'Should be \'type name\'.');
        continue;
      }

      String type = tk[0];
      String name = tk[1];

      if (type == 'float') {
        if (value is num) {
          value = [value.toDouble()];
        }
        addFloat(name, value);
      } else if (type == 'int' || type == 'integer') {
        if (value is num) {
          value = [value.toInt()];
        }
        addInt(name, value);
      } else if (type == 'bool' || type == 'boolean') {
        if (value is bool) {
          value = [value];
        }
        addBool(name, value);
      } else if (type == 'point') {
        addPoint(name, value);
      } else if (type == 'vector') {
        addVector(name, value);
      } else if (type == 'normal') {
        addNormal(name, value);
      } else if (type == 'string') {
        addString(name, value);
      } else if (type == 'string') {
        addString(name, value);
      } else if (type == 'rgb' || type == 'color') {
        addRGBSpectrum(name, value);
      } else if (type == 'xyz') {
        addXYZSpectrum(name, value);
      } else {
        LogError('Unhandled parameter type: $type');
      }
    }
  }

  Map toJson() {
    Map json = {};
    for (var p in bools) {
      json['bool ${p.name}'] = p.data;
    }

    for (var p in ints) {
      json['int ${p.name}'] = p.data;
    }

    for (var p in floats) {
      json['float ${p.name}'] = p.data;
    }

    for (var p in points) {
      json['point ${p.name}'] = p.data;
    }

    for (var p in vectors) {
      json['vector ${p.name}'] = p.data;
    }

    for (var p in normals) {
      json['normal ${p.name}'] = p.data;
    }

    for (var p in spectra) {
      json['color ${p.name}'] = p.data;
    }

    for (var p in strings) {
      json['string ${p.name}'] = p.data;
    }

    return json;
  }

  void addFloat(String name, List<double> data) {
    name = name.toLowerCase();
    eraseFloat(name);
    floats.add(new ParamSetItem<double>(name, data));
  }

  void addInt(String name, List<int> data) {
    name = name.toLowerCase();
    eraseInt(name);
    ints.add(new ParamSetItem<int>(name, data));
  }

  void addBool(String name, List<bool> data) {
    name = name.toLowerCase();
    eraseBool(name);
    bools.add(new ParamSetItem<bool>(name, data));
  }

  void addPoint(String name, List data) {
    name = name.toLowerCase();
    erasePoint(name);
    if (data[0] is num) {
      int numPts = data.length ~/ 3;
      List<Point> pts = new List<Point>(numPts);
      for (int i = 0, j = 0; i < numPts; ++i, j += 3) {
        pts[i] = new Point(data[j], data[j + 1], data[j + 2]);
      }
      points.add(new ParamSetItem<Point>(name, pts));
    } else if (data[0] is Point) {
      points.add(new ParamSetItem<Point>(name, data));
    }
  }

  void addVector(String name, List data) {
    name = name.toLowerCase();
    eraseVector(name);
    if (data[0] is num) {
      int numVecs = data.length ~/ 3;
      List<Vector> vecs = new List<Vector>(numVecs);
      for (int i = 0, j = 0; i < numVecs; ++i, j += 3) {
        vecs[i] = new Vector(data[j], data[j + 1], data[j + 2]);
      }
      vectors.add(new ParamSetItem<Vector>(name, vecs));
    } else if (data[0] is Vector) {
      vectors.add(new ParamSetItem<Vector>(name, data));
    }
  }

  void addNormal(String name, List data) {
    name = name.toLowerCase();
    eraseNormal(name);
    if (data[0] is Normal) {
      normals.add(new ParamSetItem<Normal>(name, data));
    } else if (data[0] is num) {
      int numNorms = data.length ~/ 3;
      List<Normal> norms = new List<Normal>(numNorms);
      for (int i = 0, j = 0; i < numNorms; ++i, j += 3) {
        norms[i] = new Normal(data[j], data[j + 1], data[j + 2]);
      }
      normals.add(new ParamSetItem<Normal>(name, norms));
    }
  }

  void addString(String name, List<String> data) {
    name = name.toLowerCase();
    eraseString(name);
    strings.add(new ParamSetItem<String>(name, data));
  }

  void addTexture(String name, List<String> data) {
    name = name.toLowerCase();
    eraseTexture(name);
    textures.add(new ParamSetItem<String>(name, data));
  }

  void addRGBSpectrum(String name, List<double> data) {
    name = name.toLowerCase();
    eraseSpectrum(name);
    assert(data.length % 3 == 0);
    int nItems = data.length ~/ 3;
    List<Spectrum> s = new List<Spectrum>(nItems);
    for (int i = 0, di = 0; i < nItems; ++i, di += 3) {
      s[i] = new Spectrum.rgb(data[di], data[di + 1], data[di + 2]);
    }
    spectra.add(new ParamSetItem<Spectrum>(name, s));
  }

  void addXYZSpectrum(String name, List<double> data) {
    name = name.toLowerCase();
    eraseSpectrum(name);
    assert(data.length % 3 == 0);
    int nItems = data.length ~/ 3;
    List<Spectrum> s = new List<Spectrum>(nItems);
    for (int i = 0, di = 0; i < nItems; ++i, di += 3) {
      s[i] = new Spectrum.xyz(data[di], data[di + 1], data[di + 2]);
    }
    spectra.add(new ParamSetItem<Spectrum>(name, s));
  }

  void addBlackbodySpectrum(String name, List<double> data) {
    name = name.toLowerCase();
    eraseSpectrum(name);
    assert(data.length % 2 == 0);
    int nItems = data.length ~/ 2;
    List<Spectrum> s = new List<Spectrum>(nItems);
    Float32List bb = new Float32List(Spectrum.NUM_CIE_SAMPLES);
    for (int i = 0, di = 0; i < nItems; ++i, di += 2) {
      Spectrum.Blackbody(Spectrum.CIE_lambda, data[di], bb);

      SampledSpectrum bs =
          new SampledSpectrum.fromSampled(Spectrum.CIE_lambda, bb);

      s[i] = new Spectrum.from(bs) * data[di + 1];
    }
    spectra.add(new ParamSetItem<Spectrum>(name, s));
  }

  void addSpectrumFiles(String name, List<String> filenames,
                        List<Future> futures) {
    List<Spectrum> s = new List<Spectrum>(filenames.length);

    List<Future> subFutures = [];
    for (int fi = 0; fi < filenames.length; ++fi) {
      String path = filenames[fi];

      if (ResourceManager.HasResource(path)) {
        var bytes = ResourceManager.GetResource(path);
        if (bytes is! List<int>) {
          continue;
        }

        List<double> values = ReadFloatFile(bytes, path);
        int numSamples = values.length ~/ 2;
        Float32List wls = new Float32List(numSamples);
        Float32List v = new Float32List(numSamples);
        for (int i = 0, j = 0; i < numSamples; ++i) {
          wls[i] = values[j++];
          v[i] = values[j++];
        }

        s[fi] = new Spectrum.fromSampled(wls, v);
      } else {
        if (futures == null) {
          LogDebug('UNABLE TO LOAD SPECTRUM FILE $path');
        } else {
          LogDebug('LOADING SPECTRUM FILE $path');

          Completer sc = new Completer();
          ResourceManager.RequestFile(path).then((bytes) {
            LogDebug('FINISHED SPECTRUM FILE $path');
            List<double> values = ReadFloatFile(bytes, path);
            int numSamples = values.length ~/ 2;
            Float32List wls = new Float32List(numSamples);
            Float32List v = new Float32List(numSamples);
            for (int i = 0, j = 0; i < numSamples; ++i) {
              wls[i] = values[j++];
              v[i] = values[j++];
            }
            s[fi] = new Spectrum.fromSampled(wls, v);
            sc.complete();
          });

          subFutures.add(sc.future);
        }
      }
    }

    if (futures == null) {
      spectra.add(new ParamSetItem<Spectrum>(name, s));
    } else {
      if (subFutures.isEmpty) {
        spectra.add(new ParamSetItem<Spectrum>(name, s));
      } else {
        Completer c = new Completer();
        futures.add(c.future);

        Future.wait(subFutures).then((e) {
          spectra.add(new ParamSetItem<Spectrum>(name, s));
          c.complete();
        });
      }
    }
  }

  void addSampledSpectrum(String name, List<double> data) {
    eraseSpectrum(name);
    int nItems = data.length;
    assert(nItems % 2 == 0);
    nItems ~/= 2;

    Float32List wl = new Float32List(nItems);
    Float32List v = new Float32List(nItems);

    for (int i = 0, j = 0; i < nItems; ++i) {
      wl[i] = data[j++];
      v[i] = data[j++];
    }

    List<Spectrum> s = [new Spectrum.fromSampled(wl, v)];

    spectra.add(new ParamSetItem<Spectrum>(name, s));
  }

  bool eraseInt(String n) {
    for (int i = 0; i < ints.length; ++i) {
      if (ints[i].name == n) {
        ints.removeAt(i);
        return true;
      }
    }
    return false;
  }

  bool eraseBool(String n) {
    for (int i = 0; i < bools.length; ++i) {
      if (bools[i].name == n) {
        bools.removeAt(i);
        return true;
      }
    }
    return false;
  }

  bool eraseFloat(String n) {
    for (int i = 0; i < floats.length; ++i) {
      if (floats[i].name == n) {
        floats.removeAt(i);
        return true;
      }
    }
    return false;
  }

  bool erasePoint(String n) {
    for (int i = 0; i < points.length; ++i) {
      if (points[i].name == n) {
        points.removeAt(i);
        return true;
      }
    }
    return false;
  }

  bool eraseVector(String n) {
    for (int i = 0; i < vectors.length; ++i) {
      if (vectors[i].name == n) {
        vectors.removeAt(i);
        return true;
      }
    }
    return false;
  }

  bool eraseNormal(String n) {
    for (int i = 0; i < normals.length; ++i) {
      if (normals[i].name == n) {
        normals.removeAt(i);
        return true;
      }
    }
    return false;
  }

  bool eraseSpectrum(String n) {
    for (int i = 0; i < spectra.length; ++i) {
      if (spectra[i].name == n) {
        spectra.removeAt(i);
        return true;
      }
    }
    return false;
  }

  bool eraseString(String n) {
    for (int i = 0; i < strings.length; ++i) {
      if (strings[i].name == n) {
        strings.removeAt(i);
        return true;
      }
    }
    return false;
  }

  bool eraseTexture(String n) {
    for (int i = 0; i < textures.length; ++i) {
      if (textures[i].name == n) {
        textures.removeAt(i);
        return true;
      }
    }
    return false;
  }


  double findOneFloat(String name, double d) {
    name = name.toLowerCase();
    for (int i = 0; i < floats.length; ++i) {
      if (floats[i].name == name && floats[i].data.length == 1) {
        floats[i].lookedUp = true;
        return floats[i].data[0];
      }
    }
    return d;
  }

  int findOneInt(String name, int d) {
    name = name.toLowerCase();
    for (int i = 0; i < ints.length; ++i) {
      if (ints[i].name == name && ints[i].data.length == 1) {
        ints[i].lookedUp = true;
        return ints[i].data[0];
      }
    }
    return d;
  }

  bool findOneBool(String name, bool d) {
    name = name.toLowerCase();
    for (int i = 0; i < bools.length; ++i) {
      if (bools[i].name == name && bools[i].data.length == 1) {
        bools[i].lookedUp = true;
        return bools[i].data[0];
      }
    }
    return d;
  }

  Point findOnePoint(String name, Point d) {
    name = name.toLowerCase();
    for (int i = 0; i < points.length; ++i) {
      if (points[i].name == name && points[i].data.length == 1) {
        points[i].lookedUp = true;
        return points[i].data[0];
      }
    }
    return d;
  }

  Vector findOneVector(String name, Vector d) {
    name = name.toLowerCase();
    for (int i = 0; i < vectors.length; ++i) {
      if (vectors[i].name == name && vectors[i].data.length == 1) {
        vectors[i].lookedUp = true;
        return vectors[i].data[0];
      }
    }
    return d;
  }

  Normal findOneNormal(String name, Normal d) {
    name = name.toLowerCase();
    for (int i = 0; i < normals.length; ++i) {
      if (normals[i].name == name && normals[i].data.length == 1) {
        normals[i].lookedUp = true;
        return normals[i].data[0];
      }
    }
    return d;
  }

  Spectrum findOneSpectrum(String name, Spectrum d) {
    name = name.toLowerCase();
    for (int i = 0; i < spectra.length; ++i) {
      if (spectra[i].name == name) {
        spectra[i].lookedUp = true;

        if (spectra[i].data == null) {
          return d;
        }

        if (spectra[i].data.length == 1) {
          spectra[i].lookedUp = true;
          return spectra[i].data[0];
        }
      }
    }
    return d;
  }

  String findOneString(String name, String d) {
    name = name.toLowerCase();
    for (int i = 0; i < strings.length; ++i) {
      if (strings[i].name == name && strings[i].data.length == 1) {
        strings[i].lookedUp = true;
        return strings[i].data[0];
      }
    }
    return d;
  }

  String findOneFilename(String name, String d) {
    name = name.toLowerCase();
    String filename = findOneString(name, "");
    if (filename == "") {
      return d;
    }
    return filename;
  }

  String findTexture(String name) {
    name = name.toLowerCase();
    String d = "";
    for (int i = 0; i < textures.length; ++i) {
      if (textures[i].name == name &&
          textures[i].data.length == 1) {
        textures[i].lookedUp = true;
        return textures[i].data[0];
      }
    }
    return d;
  }

  List<double> findFloat(String name) {
    name = name.toLowerCase();
    for (int i = 0; i < floats.length; ++i) {
      if (floats[i].name == name) {
        floats[i].lookedUp = true;
        return floats[i].data;
      }
    }
    return null;
  }

  List<int> findInt(String name) {
    name = name.toLowerCase();
    for (int i = 0; i < ints.length; ++i) {
      if (ints[i].name == name) {
        ints[i].lookedUp = true;
        return ints[i].data;
      }
    }
    return null;
  }

  List<bool> findBool(String name) {
    name = name.toLowerCase();
    for (int i = 0; i < bools.length; ++i) {
      if (bools[i].name == name) {
        bools[i].lookedUp = true;
        return bools[i].data;
      }
    }
    return null;
  }

  List<Point> findPoint(String name) {
    name = name.toLowerCase();
    for (int i = 0; i < points.length; ++i) {
      if (points[i].name == name) {
        points[i].lookedUp = true;
        return points[i].data;
      }
    }
    return null;
  }

  List<Vector> findVector(String name) {
    name = name.toLowerCase();
    for (int i = 0; i < vectors.length; ++i) {
      if (vectors[i].name == name) {
        vectors[i].lookedUp = true;
        return vectors[i].data;
      }
    }
    return null;
  }

  List<Normal> findNormal(String name) {
    name = name.toLowerCase();
    for (int i = 0; i < normals.length; ++i) {
      if (normals[i].name == name) {
        normals[i].lookedUp = true;
        return normals[i].data;
      }
    }
    return null;
  }

  List<Spectrum> findSpectrum(String name) {
    name = name.toLowerCase();
    for (int i = 0; i < spectra.length; ++i) {
      if (spectra[i].name == name) {
        spectra[i].lookedUp = true;
        return spectra[i].data;
      }
    }
    return null;
  }

  List<String> findString(String name) {
    name = name.toLowerCase();
    for (int i = 0; i < strings.length; ++i) {
      if (strings[i].name == name) {
        strings[i].lookedUp = true;
        return strings[i].data;
      }
    }
    return null;
  }

  void reportUnused() {
    _reportUnused(bools);
    _reportUnused(ints);
    _reportUnused(floats);
    _reportUnused(points);
    _reportUnused(vectors);
    _reportUnused(normals);
    _reportUnused(spectra);
    _reportUnused(strings);
    _reportUnused(textures);
  }

  void _reportUnused(List<ParamSetItem> list) {
    for (int i = 0; i < list.length; ++i) {
      if (!list[i].lookedUp) {
        LogWarning('Parameter ${list[i].name} not used');
      }
    }
  }

  void clear() {
    bools.clear();
    ints.clear();
    floats.clear();
    points.clear();
    vectors.clear();
    normals.clear();
    spectra.clear();
    strings.clear();
    textures.clear();
  }

  String toString() {
    String out = '';
    for (var p in bools) {
      out += '"bool ${p.name}" [';
      for (int i = 0; i < p.data.length; ++i) {
        if (i != 0) {
          out += ' ';
        }
        out += '${p.data[i] ? 1 : 0}';
      }
      out += '] ';
    }

    for (var p in ints) {
      out += _paramToString('integer', p);
    }

    for (var p in floats) {
      out += _paramToString('float', p);
    }

    for (var p in points) {
      out += _paramToString('point', p);
    }

    for (var p in vectors) {
      out += _paramToString('vector', p);
    }

    for (var p in normals) {
      out += _paramToString('normals', p);
    }

    for (var p in spectra) {
      out += _paramToString('color', p);
    }

    for (var p in strings) {
      out += '"string ${p.name}" [';
      for (int i = 0; i < p.data.length; ++i) {
        if (i != 0) {
          out += ' ';
        }
        out += '"${p.data[i]}"';
      }
      out += '] ';
    }

    for (var p in textures) {
      out += '"texture ${p.name}" [';
      for (int i = 0; i < p.data.length; ++i) {
        if (i != 0) {
          out += ' ';
        }
        out += '"${p.data[i]}"';
      }
      out += '] ';
    }

    return out;
  }

  String _paramToString(String type, ParamSetItem item) {
    String out;
    out = '"$type ${item.name}" [';
    for (int i = 0; i < item.data.length; ++i) {
      if (i != 0) {
        out += ' ';
      }
      out += '${item.data[i]}';
    }
    out += '] ';
    return out;
  }

  List<ParamSetItem<bool>> bools = [];
  List<ParamSetItem<int>> ints = [];
  List<ParamSetItem<double>> floats = [];
  List<ParamSetItem<Point>> points = [];
  List<ParamSetItem<Vector>> vectors = [];
  List<ParamSetItem<Normal>> normals = [];
  List<ParamSetItem<Spectrum>> spectra = [];
  List<ParamSetItem<String>> strings = [];
  List<ParamSetItem<String>> textures = [];
}


class ParamSetItem<T> {
  ParamSetItem(this.name, this.data);

  String name;
  List<T> data;
  bool lookedUp = false;
}
