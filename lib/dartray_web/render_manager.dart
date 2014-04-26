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
part of dartray_web;

/**
 * [RenderManager] provides the interface for loading scenes and rendering.
 * This implementation is used by dart:http web apps, and uses HttpRequest
 * calls to load resources such as geometry and texture images.
 */
class RenderManager extends RenderManagerInterface {
  RenderManager([String scenePath = '']) :
    super(scenePath);

  Future<List<int>> loadFile(String path, [Future future]) {
    if (future != null) {
      futures.add(future);
    }

    Completer<List<int>> c = new Completer<List<int>>();

    _loadFile(path).then((bytes) {
      c.complete(bytes);
    }).catchError((e) {
      LogError(e.toString());
      c.complete(null);
    });

    return c.future;
  }

  Future<List<int>> _loadFile(String path) {
    Completer<List<int>> c = new Completer<List<int>>();

    // As of dart 1.3.0, we can't call HttpRequest from an isolate. To work
    // around this, we send a request message back to the main isolate asking
    // it to load the file, and it will send the file data back here.
    if (isolate != null) {
      isolate.requestResponse({'cmd': 'file', 'path': path}).then((bytes) {
        if (bytes is String) {
          c.complete(bytes.codeUnits);
          return;
        } else if (bytes is ByteBuffer) {
          c.complete(new Uint8List.view(bytes));
          return;
        } else if (bytes is List<int>) {
          c.complete(bytes);
          return;
        } else {
          LogError('Unknown HttpRequest response type');
        }
      });
    } else {
      path = scenePath + '/' + path;

      Html.HttpRequest.request(path, method: 'GET',
                               mimeType: 'text\/plain; charset=x-user-defined')
      .then((resp) {
        if (resp.response is String) {
          String s = resp.response;
          Uint8List bytes = new Uint8List.fromList(s.codeUnits);
          c.complete(bytes);
          return;
        } else if (resp.response is ByteBuffer) {
          c.complete(new Uint8List.view(resp.response));
          return;
        } else if (resp.response is List<int>) {
          c.complete(resp.response);
          return;
        } else {
          LogError('Unknown HttpRequest response type');
        }
      }).catchError((e) {
        LogError('Error Loading Resource $path');
      });
    }

    return c.future;
  }
}
