/****************************************************************************
 *  Copyright (C) 2014 by Brendan Duncan.                                   *
 *                                                                          *
 *  This file is part of DartRay.                                           *
 *                                                                          *
 *  Licensed under the Apache License, Version 2.0 (the "License");         *
 *  you may not use this file except in compliance with the License.        *
 *  You may obtain a copy of the License at                                 *
 *                                                                          *
 *  http://www.apache.org/licenses/LICENSE-2.0                              *
 *                                                                          *
 *  Unless required by applicable law or agreed to in writing, software     *
 *  distributed under the License is distributed on an "AS IS" BASIS,       *
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 *  See the License for the specific language governing permissions and     *
 *  limitations under the License.                                          *
 *                                                                          *
 *   This project is based on PBRT v2 ; see http://www.pbrt.org             *
 *   pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.*
 ****************************************************************************/
part of pbrt;

class PbrtParser {
  Pbrt pbrt;

  PbrtParser(this.pbrt);

  Future parse(String input) {
    Completer c = new Completer();
    _loadIncludes(input).then((x) {
      _parse(input);
      c.complete();
    });

    return c.future;
  }

  Future _loadIncludes(String input) {
    PbrtLexer _lexer = new PbrtLexer(input);

    List<Future<List<int>>> futures = [];

    Completer c = new Completer();

    int tk = _lexer.nextToken();
    while (!_lexer.isEof()) {
      Map cmd = _parseCommand(_lexer);
      if (cmd == null) {
        _lexer.nextToken();
        continue;
      }

      String name = cmd['name'].toLowerCase();

      switch (name) {
        case 'include':
          futures.add(pbrt.resourceManager.requrestFile(cmd['value']));
          break;
      }
    }

    Future.wait(futures)
          .then((List responses) {
            c.complete();
          })
          .catchError((e) {
            LogError(e);
            c.completeError(e);
          });

    return c.future;
  }

  void _parse(String input) {
    PbrtLexer _lexer = new PbrtLexer(input);

    int tk = _lexer.nextToken();
    while (!_lexer.isEof()) {
      Map cmd = _parseCommand(_lexer);

      if (cmd == null) {
        _lexer.nextToken();
        continue;
      }

      String name = cmd['name'].toLowerCase();

      switch (name) {
        case 'accelerator':
          pbrt.accelerator(cmd['type'], cmd['params']);
          break;
        case 'activetransform':
          if (cmd['type'] == 'StartTime') {
            pbrt.activeTransformStartTime();
          } else if (cmd['type'] == 'EndTime') {
            pbrt.activeTransformEndTime();
          } else if (cmd['type'] == 'All') {
            pbrt.activeTransformAll();
          }
          break;
        case 'arealightsource':
          pbrt.areaLightSource(cmd['type'], cmd['params']);
          break;
        case 'attributebegin':
          pbrt.attributeBegin();
          break;
        case 'attributeend':
          pbrt.attributeEnd();
          break;
        case 'camera':
          pbrt.camera(cmd['type'], cmd['params']);
          break;
        case 'coordinatesystem':
          pbrt.coordinateSystem(cmd['type']);
          break;
        case 'coordsystransform':
          pbrt.coordSysTransform(cmd['type']);
          break;
        case 'concattransform':
          if (!cmd.containsKey('values') || cmd['values'].length != 16) {
            LogWarning('ConcatTransform requires 16 values, ${cmd['values'].length} found.');
          } else {
            var v = cmd['values'];
            Matrix4x4 m = new Matrix4x4.fromList(v);
            pbrt.concatTransform(m);
          }
          break;
        case 'film':
          pbrt.film(cmd['type'], cmd['params']);
          break;
        case 'ident':
          pbrt.identity();
          break;
        case 'lightsource':
          pbrt.lightSource(cmd['type'], cmd['params']);
          break;
        case 'lookat':
          if (!cmd.containsKey('values') || cmd['values'].length != 9) {
            LogWarning('LookAt requires 9 values');
          } else {
            var v = cmd['values'];
            pbrt.lookAt(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
          }
          break;
        case 'makenamedmaterial':
          pbrt.makeNamedMaterial(cmd['id'], cmd['params']);
          break;
        case 'material':
          pbrt.material(cmd['type'], cmd['params']);
          break;
        case 'namedmaterial':
          pbrt.namedMaterial(cmd['type']);
          break;
        case 'objectbegin':
          pbrt.objectBegin(cmd['type']);
          break;
        case 'objectend':
          pbrt.objectEnd();
          break;
        case 'objectinstance':
          pbrt.objectInstance(cmd['type']);
          break;
        case 'pixelfilter':
          pbrt.pixelFilter(cmd['type'], cmd['params']);
          break;
        case 'renderer':
          pbrt.renderer(cmd['type'], cmd['params']);
          break;
        case 'reverseorientation':
          pbrt.reverseOrientation();
          break;
        case 'rotate':
          if (!cmd.containsKey('values') || cmd['values'].length != 4) {
            LogWarning('Rotate requires 4 values');
          } else {
            var v = cmd['values'];
            pbrt.rotate(v[0], v[1], v[2], v[3]);
          }
          break;
        case 'sampler':
          pbrt.sampler(cmd['type'], cmd['params']);
          break;
        case 'scale':
          if (!cmd.containsKey('values') || cmd['values'].length != 3) {
            LogWarning('Scale requires 3 values');
          } else {
            var v = cmd['values'];
            pbrt.scale(v[0], v[1], v[2]);
          }
          break;
        case 'shape':
          pbrt.shape(cmd['type'], cmd['params']);
          break;
        case 'surfaceintegrator':
          pbrt.surfaceIntegrator(cmd['type'], cmd['params']);
          break;
        case 'texture':
          pbrt.texture(cmd['id'], cmd['type'], cmd['class'], cmd['params']);
          break;
        case 'translate':
          if (!cmd.containsKey('values') || cmd['values'].length != 3) {
            LogWarning('Translate requires 3 values');
          } else {
            var v = cmd['values'];
            pbrt.translate(v[0], v[1], v[2]);
          }
          break;
        case 'transform':
          if (!cmd.containsKey('values') || cmd['values'].length != 16) {
            LogWarning('Transform requires 16 values');
          } else {
            var v = cmd['values'];
            Matrix4x4 m = new Matrix4x4.fromList(v);
            pbrt.transform(m);
          }
          break;
        case 'transformbegin':
          pbrt.transformBegin();
          break;
        case 'transformend':
          pbrt.transformEnd();
          break;
        case 'transformtimes':
          if (!cmd.containsKey('values') || cmd['values'].length != 2) {
            LogWarning('TransformTimes requires 2 values');
          } else {
            var v = cmd['values'];
            pbrt.transformTimes(v[0], v[1]);
          }
          break;
        case 'volume':
          pbrt.volume(cmd['type'], cmd['params']);
          break;
        case 'volumeintegrator':
          pbrt.volumeIntegrator(cmd['type'], cmd['params']);
          break;
        case 'worldbegin':
          pbrt.worldBegin();
          break;
        case 'worldend':
          pbrt.worldEnd();
          break;
        case 'include':
          List<int> bytes = pbrt.resourceManager.getFile(cmd['value']);
          if (bytes != null) {
            String inc = new String.fromCharCodes(bytes);
            _parse(inc);
          }
          break;
        default:
          LogWarning('Unhandled command ${cmd['name']}');
          break;
      }
    }
  }

  Map _parseCommand(PbrtLexer _lexer) {
    if (_lexer.currentToken != PbrtLexer.TOKEN_IDENTIFIER) {
      return null;
    }

    var cmd = {};
    cmd['name'] = _lexer.currentTokenString;

    String name = cmd['name'].toLowerCase();

    if (name == 'include') {
      int tk = _lexer.nextToken();
      cmd['value'] = _lexer.currentTokenString;
      return cmd;
    }

    if (name == 'activetransform') {
      int tk = _lexer.nextToken();
      cmd['type'] = _lexer.currentTokenString;
      return cmd;
    }

    List params = [];

    int tk;
    if (name == 'texture') {
      tk = _lexer.nextToken();
      cmd['id'] = _lexer.currentTokenString;

      tk = _lexer.nextToken();
      cmd['type'] = _lexer.currentTokenString;

      tk = _lexer.nextToken();
      cmd['class'] = _lexer.currentTokenString;

      tk = _lexer.nextToken();
    } else if (name == 'makenamedmaterial') {
      tk = _lexer.nextToken();
      cmd['id'] = _lexer.currentTokenString;

      tk = _lexer.nextToken();
      params.add({'type': 'string',
                  'name': 'type',
                  'value': [_lexer.currentTokenString]});

      tk = _lexer.nextToken();
    } else {
      tk = _lexer.nextToken();
      if (tk == PbrtLexer.TOKEN_IDENTIFIER) {
        return cmd;
      }

      if (tk == PbrtLexer.TOKEN_STRING) {
        cmd['type'] = _lexer.currentTokenString;
        tk = _lexer.nextToken();
      } else if (tk == PbrtLexer.TOKEN_LEFT_BRACKET) {
        List values = [];
        tk = _lexer.nextToken();
        while (tk != PbrtLexer.TOKEN_RIGHT_BRACKET && !_lexer.isEof()) {
          if (_lexer.value != null) {
            values.add(_lexer.value);
          } else {
            values.add(_lexer.currentTokenString);
          }
          tk = _lexer.nextToken();
        }
        cmd['values'] = values;
      } else if (tk == PbrtLexer.TOKEN_NUMBER) {
        List values = [_lexer.value];
        tk = _lexer.nextToken();
        while (tk == PbrtLexer.TOKEN_NUMBER && !_lexer.isEof()) {
          values.add(_lexer.value);
          tk = _lexer.nextToken();
        }
        cmd['values'] = values;
      }
    }

    while (tk == PbrtLexer.TOKEN_STRING && !_lexer.isEof()) {
      String paramTypeName = _lexer.currentTokenString;
      tk = _lexer.nextToken();
      var paramValue;
      if (tk == PbrtLexer.TOKEN_LEFT_BRACKET) {
        List array = [];
        tk = _lexer.nextToken();
        while (tk != PbrtLexer.TOKEN_RIGHT_BRACKET && !_lexer.isEof()) {
          array.add(_lexer.currentTokenString);
          tk = _lexer.nextToken();
        }
        paramValue = array;
      } else {
        paramValue = [_lexer.currentTokenString];
      }

      List<String> paramTk = paramTypeName.split(' ');
      paramTk.removeWhere((e) => e.isEmpty);
      if (paramTk.length != 2) {
        LogWarning('Expected Parameter "Type Name", found: "$paramTypeName"');
        return cmd;
      }

      params.add({'type': paramTk[0],
                  'name': paramTk[1],
                  'value': paramValue});

      tk = _lexer.nextToken();
    }

    if (params.isNotEmpty) {
      cmd['params'] = _parseParameters(params);
    } else {
      cmd['params'] = new ParamSet();
    }

    return cmd;
  }

  ParamSet _parseParameters(List params) {
    ParamSet ps = new ParamSet();

    for (var p in params) {
      switch (p['type']) {
        case 'float':
          var v = p['value'];
          for (int i = 0, l = v.length; i < l; ++i) {
            v[i] = double.parse(v[i]);
          }
          ps.addFloat(p['name'], v);
          break;
        case 'integer':
          var v = p['value'];
          for (int i = 0, l = v.length; i < l; ++i) {
            v[i] = int.parse(v[i]);
          }
          ps.addInt(p['name'], v);
          break;
        case 'color':
          var v = p['value'];
          for (int i = 0, l = v.length; i < l; ++i) {
            v[i] = double.parse(v[i]);
          }
          ps.addRGBSpectrum(p['name'], v);
          break;
        case 'point':
          var v = p['value'];
          for (int i = 0, l = v.length; i < l; ++i) {
            v[i] = double.parse(v[i]);
          }
          ps.addPoint(p['name'], v);
          break;
        case 'normal':
          var v = p['value'];
          for (int i = 0, l = v.length; i < l; ++i) {
            v[i] = double.parse(v[i]);
          }
          ps.addNormal(p['name'], v);
          break;
        case 'vector':
          var v = p['value'];
          for (int i = 0, l = v.length; i < l; ++i) {
            v[i] = double.parse(v[i]);
          }
          ps.addVector(p['name'], v);
          break;
        case 'string':
          ps.addString(p['name'], p['value']);
          break;
        case 'texture':
          ps.addTexture(p['name'], p['value']);
          break;
        case 'spectrum':
          var v = p['value'];
          if (v.isNotEmpty) {
            if (v[0] is String) {
              ps.addSpectrumFiles(p['name'], v);
            } else {
              if ((v.length % 2) != 0) {
                LogWarning("Non-even number of values given with sampled spectrum "
                           "parameter \"${p['name']}\". Ignoring extra.");
              }
              for (int i = 0, l = v.length; i < l; ++i) {
                v[i] = double.parse(v[i]);
              }
              ps.addSampledSpectrum(p['name'], v);
            }
          }
          break;
        case 'blackbody':
          var v = p['value'];
          if (v.isNotEmpty) {
            for (int i = 0, l = v.length; i < l; ++i) {
              v[i] = double.parse(v[i]);
            }
            ps.addBlackbodySpectrum(p['name'], v);
          }
          break;
        case 'bool':
          var v = p['value'];
          for (int i = 0, l = v.length; i < l; ++i) {
            if (v[i].toLowerCase() == 'false') {
              v[i] = false;
            } else if (v[i].toLowerCase() == 'true') {
              v[i] = true;
            } else {
              LogWarning('Invalid boolean value "${v[i]}".');
              v[i] = false;
            }
          }
          ps.addBool(p['name'], v);
          break;
        default:
          LogWarning('Unhandled parameter type ${p['type']}');
          break;
      }
    }

    return ps;
  }
}
