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
part of dartray_io;

class RenderManager extends RenderManagerInterface {
  RenderManager([String scenePath = '']) :
    super(scenePath);

  Future<List<int>> loadFile(String file) {
    Completer<List<int>> c = new Completer<List<int>>();

    String path = _findFile(file);

    if (path == null) {
      LogInfo('File not found: $file');
      c.complete(null);
      return c.future;
    }

    LogDebug('LOAD $path');
    Future<List<int>> f = new File(path).readAsBytes();
    f.then((bytes) {
      c.complete(bytes);
    });

    return c.future;
  }

  String _findFile(String file) {
    if (FileSystemEntity.isFileSync(file)) {
      return file;
    }

    if (FileSystemEntity.isFileSync(scenePath + '/' + file)) {
      return scenePath + '/' + file;
    }

    for (String p in includePaths) {
      if (FileSystemEntity.isFileSync(p + '/' + file)) {
        return p + '/' + file;
      }
    }

    return null;
  }

  List<String> includePaths = [];
}
