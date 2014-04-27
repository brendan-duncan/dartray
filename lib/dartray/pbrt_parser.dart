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
part of dartray;

class PbrtParser {
  DartRay dartray;
  ResourceManager resourceManager;

  PbrtParser(DartRay pbrt) :
    this.dartray = pbrt,
    resourceManager = pbrt.resourceManager;

  Future parse(String file) {
    Stopwatch t = new Stopwatch()..start();

    LogInfo('LOADING Scene $file');
    Completer c = new Completer();

    resourceManager.requestFile(file).then((List<int> input) {
      input = _decodeFile(file, input);
      LogDebug('FINISHED Loading $file. Scanning for includes.');
      _loadIncludes(input, file).then((_) {
        LogDebug('FINISHED Includes.');
        _parse(input, file).then((_) {
          t.stop();
          LogInfo('FINISHED Parsing Scene: ${t.elapsed}');
          c.complete();
        });
      });
    });

    return c.future;
  }

  bool _isSceneFileName(String file) {
    int i = file.lastIndexOf('/');
    if (i == -1) {
      i = file.lastIndexOf('\\');
    }

    if (i != -1) {
      file = file.substring(i + 1);
    }

    return file == 'scene.pbrt' || file == 'scene.pbrt.z' ||
           file == 'scene.pbrt.gz' || file == 'scene.pbrt.bz2';
  }

  bool _isArchive(String file) {
    return file.endsWith('.zip') || file.endsWith('.tar') ||
           file.endsWith('.tgz') || file.endsWith('.tbz') ||
           file.endsWith('.tar.gz') || file.endsWith('.tar.bz2');
  }

  Future _loadIncludes(List<int> input, String path) {
    PbrtLexer _lexer = new PbrtLexer(input, path);

    List<Future> futures = [];
    List<String> paths = [];

    Completer c = new Completer();
    int tk = _lexer.nextToken();
    while (!_lexer.isEof()) {
      Map cmd = _parseCommand(_lexer, futures);
      if (cmd == null) {
        _lexer.nextToken();
        continue;
      }

      String name = cmd['name'].toLowerCase();

      if (name == 'include') {
        String file = cmd['value'];
        if (!dartray.resourceManager.hasResource(file)) {
          LogDebug('LOADING INCLUDE ${file}');
          Future f = dartray.resourceManager.requestFile(file);
          futures.add(f);
          paths.add(cmd['value']);
          f.then((List<int> input) {
            _decodeFile(file, input);
            LogDebug('FINISHED LOADING INCLUDE ${cmd['value']}');
          });
        }
      }
    }

    if (futures.isNotEmpty) {
      Future.wait(futures).then((List responses) {
        LogDebug('LOADING SUB-INCLUDES');
        List<Future<List<int>>> subFutures = [];

        if (responses.isNotEmpty) {
          for (int i = 0; i < responses.length; ++i) {
            List<int> inc = responses[i];
            if (inc != null && inc.isNotEmpty) {
              subFutures.add(_loadIncludes(inc, '@'));
            }
          }

          Future.wait(subFutures).then((List subResponses) {
            LogDebug('FINISHED SUB-INCLUDES');
            c.complete();
          }).catchError((e) {
            LogError(e);
            c.completeError(e);
          });
        } else {
          c.complete();
        }
      }).catchError((e) {
        LogError(e.toString());
        c.completeError(e);
      });
    } else {
      c.complete();
    }

    return c.future;
  }

  Future _parse(List<int> input, String path) {
    LogDebug('PARSING FILE $path');
    Completer c = new Completer();

    PbrtLexer _lexer = new PbrtLexer(input, path);

    try {
      FutureWhileLoop(() => _nextCommand(_lexer)).then((_) {
        c.complete();
      }).catchError((e) {
        c.completeError(e);
        LogError('EXCEPTION: $e');
      });
    } catch (e) {
      c.completeError(e);
      LogError('EXCEPTION: $e');
    }

    return c.future;
  }

  Future _nextCommand(PbrtLexer lexer) {
    if (lexer.isEof()) {
      // Returning null will break the parsing loop
      return null;
    }

    Map cmd = _parseCommand(lexer, null);

    if (cmd == null) {
      // Nothing was parsed, then there's nothing to do.
      lexer.nextToken();
      Completer completer = new Completer();
      completer.complete();
      return completer.future;
    }

    String name = cmd['name'].toLowerCase();

    Future cmdFuture;

    try {
      switch (name) {
        case 'accelerator':
          dartray.accelerator(cmd['type'], cmd['params']);
          break;
        case 'activetransform':
          if (cmd['type'] == 'StartTime') {
            dartray.activeTransformStartTime();
          } else if (cmd['type'] == 'EndTime') {
            dartray.activeTransformEndTime();
          } else if (cmd['type'] == 'All') {
            dartray.activeTransformAll();
          }
          break;
        case 'arealightsource':
          dartray.areaLightSource(cmd['type'], cmd['params']);
          break;
        case 'attributebegin':
          dartray.attributeBegin();
          break;
        case 'attributeend':
          dartray.attributeEnd();
          break;
        case 'camera':
          dartray.camera(cmd['type'], cmd['params']);
          break;
        case 'coordinatesystem':
          dartray.coordinateSystem(cmd['type']);
          break;
        case 'coordsystransform':
          dartray.coordSysTransform(cmd['type']);
          break;
        case 'concattransform':
          if (!cmd.containsKey('values') || cmd['values'].length != 16) {
            LogWarning('ConcatTransform requires 16 values, '
                       '${cmd['values'].length} found.');
          } else {
            var v = cmd['values'];
            Matrix4x4 m = Matrix4x4.Transpose(new Matrix4x4.fromList(v));
            dartray.concatTransform(m);
          }
          break;
        case 'film':
          dartray.film(cmd['type'], cmd['params']);
          break;
        case 'ident':
        case 'identity':
          dartray.identity();
          break;
        case 'lightsource':
          dartray.lightSource(cmd['type'], cmd['params']);
          break;
        case 'lookat':
          if (!cmd.containsKey('values') || cmd['values'].length != 9) {
            LogWarning('LookAt requires 9 values');
          } else {
            var v = cmd['values'];
            dartray.lookAt(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
          }
          break;
        case 'makenamedmaterial':
          dartray.makeNamedMaterial(cmd['id'], cmd['params']);
          break;
        case 'material':
          dartray.material(cmd['type'], cmd['params']);
          break;
        case 'namedmaterial':
          dartray.namedMaterial(cmd['type']);
          break;
        case 'objectbegin':
          dartray.objectBegin(cmd['type']);
          break;
        case 'objectend':
          dartray.objectEnd();
          break;
        case 'objectinstance':
          dartray.objectInstance(cmd['type']);
          break;
        case 'pixelfilter':
          dartray.pixelFilter(cmd['type'], cmd['params']);
          break;
        case 'renderer':
          dartray.renderer(cmd['type'], cmd['params']);
          break;
        case 'reverseorientation':
          dartray.reverseOrientation();
          break;
        case 'rotate':
          if (!cmd.containsKey('values') || cmd['values'].length != 4) {
            LogWarning('Rotate requires 4 values');
          } else {
            var v = cmd['values'];
            dartray.rotate(v[0], v[1], v[2], v[3]);
          }
          break;
        case 'pixels':
          dartray.pixels(cmd['type'], cmd['params']);
          break;
        case 'sampler':
          dartray.sampler(cmd['type'], cmd['params']);
          break;
        case 'scale':
          if (!cmd.containsKey('values') || cmd['values'].length != 3) {
            LogWarning('Scale requires 3 values');
          } else {
            var v = cmd['values'];
            dartray.scale(v[0], v[1], v[2]);
          }
          break;
        case 'shape':
          dartray.shape(cmd['type'], cmd['params']);
          break;
        case 'surfaceintegrator':
          dartray.surfaceIntegrator(cmd['type'], cmd['params']);
          break;
        case 'texture':
          dartray.texture(cmd['id'], cmd['type'], cmd['class'], cmd['params']);
          break;
        case 'translate':
          if (!cmd.containsKey('values') || cmd['values'].length != 3) {
            LogWarning('Translate requires 3 values');
          } else {
            var v = cmd['values'];
            dartray.translate(v[0], v[1], v[2]);
          }
          break;
        case 'transform':
          if (!cmd.containsKey('values') || cmd['values'].length != 16) {
            LogWarning('Transform requires 16 values');
          } else {
            var v = cmd['values'];
            Matrix4x4 m = new Matrix4x4.values(v[0], v[4], v[8], v[12],
                                               v[1], v[5], v[9], v[13],
                                               v[2], v[6], v[10], v[14],
                                               v[3], v[7], v[11], v[15]);
            dartray.transform(m);
          }
          break;
        case 'transformbegin':
          dartray.transformBegin();
          break;
        case 'transformend':
          dartray.transformEnd();
          break;
        case 'transformtimes':
          if (!cmd.containsKey('values') || cmd['values'].length != 2) {
            LogWarning('TransformTimes requires 2 values');
          } else {
            var v = cmd['values'];
            dartray.transformTimes(v[0], v[1]);
          }
          break;
        case 'volume':
          dartray.volume(cmd['type'], cmd['params']);
          break;
        case 'volumeintegrator':
          dartray.volumeIntegrator(cmd['type'], cmd['params']);
          break;
        case 'worldbegin':
          dartray.worldBegin();
          break;
        case 'worldend':
          cmdFuture = dartray.worldEnd();
          break;
        case 'include':
          String name = cmd['value'];
          LogInfo('INCLUDE $name');
          // don't parse archive includes (their contents were added to the
          // resources during the preprocess).
          if (!_isArchive(name)) {
            var inc = dartray.resourceManager.getResource(cmd['value']);
            if (inc is List<int>) {
              lexer.addInclude(inc, cmd['value']);
            } else {
              LogWarning('MISSING include: ${cmd['value']}');
            }
          }
          break;
        default:
          LogWarning('UNHANDLED command ${cmd['name']}');
          break;
      }
    } catch (e) {
      LogError('EXCEPTION: $e');
    }

    if (cmdFuture == null) {
      Completer completer = new Completer();
      completer.complete();
      cmdFuture = completer.future;
    }

    return cmdFuture;
  }

  Map _parseCommand(PbrtLexer _lexer, List<Future> futures) {
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
      _lexer.nextToken();
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
        LogWarning('${_lexer.path} [${_lexer.line}]: '
                   'Expected Parameter "Type Name", found: '
                   '"$paramTypeName"');
        return cmd;
      }

      params.add({'type': paramTk[0],
                  'name': paramTk[1],
                  'value': paramValue});

      tk = _lexer.nextToken();
    }

    if (params.isNotEmpty) {
      cmd['params'] = _parseParameters(params, futures);
    } else {
      cmd['params'] = new ParamSet();
    }

    return cmd;
  }

  ParamSet _parseParameters(List params, List<Future> futures) {
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
        case 'rgb':
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
               ps.addSpectrumFiles(p['name'], v, futures);
            } else {
              if ((v.length % 2) != 0) {
                LogWarning('Non-even number of values given with sampled '
                           'spectrum parameter \"${p['name']}\". '
                           'Ignoring extra.');
              }
              for (int i = 0, l = v.length; i < l; ++i) {
                v[i] = double.parse(v[i]);
              }
              ps.addSampledSpectrum(p['name'], v);
            }
          }
          break;
        case 'xyz':
          var v = p['value'];
          for (int i = 0, l = v.length; i < l; ++i) {
            v[i] = double.parse(v[i]);
          }
          ps.addXYZSpectrum(p['name'], v);
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

  /**
   * Checks to see if the file is an archive (zip, tar), and extracts its
   * contents if it is.
   */
  List<int> _decodeFile(String file, List<int> input) {
    if (!_isArchive(file)) {
      return input;
    }

    try {
      if (file.endsWith('.tgz') || file.endsWith('.tar.gz')) {
        input = new GZipDecoder().decodeBytes(input);
        file += '.tar';
      } else if (file.endsWith('.tbz') || file.endsWith('.tar.bz2')) {
        input = new BZip2Decoder().decodeBytes(input);
        file += '.tar';
      }

      if (file.endsWith('.zip')) {
        LogInfo('Decoding Zip Archive: $file');
        Archive arc = new ZipDecoder().decodeBytes(input);
        input = null;
        for (ArchiveFile f in arc) {
          List<int> fc = _decodeFile(f.name, f.content);
          if (_isSceneFileName(f.name)) {
            input = f.content;
          }
          resourceManager.setResource(f.name, f.content);
        }
      } else if (file.endsWith('.tar')) {
        LogInfo('Decoding Tar Archive: $file');
        Archive arc = new TarDecoder().decodeBytes(input);
        input = null;
        for (ArchiveFile f in arc) {
          if (_isSceneFileName(f.name)) {
            input = f.content;
          }
          resourceManager.setResource(f.name, f.content);
        }
      }
    } catch (e) {
      LogError('EXCEPTION: $e');
    }

    return input;
  }
}
