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
  RenderIsolate isolate;
  List<RenderTask> isolates;

  RenderManager([String scenePath = ''])
      : super(scenePath);

  /**
   * This is called from an Isolate to initialize the RenderManager for the
   * Isolate.
   */
  void startIsolate() {
    isolate = new RenderIsolate(this);
    isolate.start();
  }

  void pause() {
    if (isolates == null) {
      return;
    }

    for (RenderTask task in isolates) {
      task.pause();
    }
  }

  void resume() {
    if (isolates == null) {
      return;
    }

    for (RenderTask task in isolates) {
      task.resume;
    }
  }

  /**
   * [isolateUri] is the path to the script that initiates the
   * RenderIsolate.  This script usually would do nothing more than:
   * import 'dart:isolate';
   * import 'package:dartray/dartray.dart';
   * void main(List<String> args, SendPort port) {
   *   new RenderIsolate().start(port);
   * }
   */
  Future<OutputImage> render(String path, {String isolate,
    LogCallback log, PreviewCallback preview,
    RenderOverrides overrides, int numThreads: 1}) {
    if (log != null) {
      Log = log;
    }

    if (path.contains('/')) {
      int i = path.lastIndexOf('/');
      scenePath = path.substring(0, i);
      path = path.substring(i + 1);
    }

    Completer<OutputImage> completer = new Completer<OutputImage>();

    if (isolate == null) {
      LogInfo('STARTING DartRay Render');
      dartray.renderScene(path, overrides: overrides).then((output) {
        completer.complete(output);
      });

      return completer.future;
    }

    int tasksRemaining = numThreads;
    isolates = new List<RenderTask>(numThreads);

    for (int i = 0; i < numThreads; ++i) {
      isolates[i] = new RenderTask(preview, i, numThreads);
      isolates[i].render(path, isolate, overrides: overrides).then((output) {
        if (numThreads > 1) {
          if (renderOutput == null ||
              renderOutput.imageWidth != output.imageWidth ||
              renderOutput.imageHeight != output.imageHeight) {
            renderOutput = new OutputImage(0, 0, output.imageWidth,
                output.imageHeight);
          }

          for (int y = 0; y < output.height; ++y) {
            int pi = ((output.yOffset + y) * renderOutput.imageWidth * 3) +
                (output.xOffset * 3);
            for (int x = 0; x < output.width; ++x, pi += 3) {
              renderOutput.rgb[pi] = output.rgb[pi];
              renderOutput.rgb[pi + 1] = output.rgb[pi + 1];
              renderOutput.rgb[pi + 2] = output.rgb[pi + 2];
            }
          }
        } else {
          renderOutput = output;
        }

        tasksRemaining--;
        if (tasksRemaining == 0) {
          completer.complete(renderOutput);
        }
      }, onError: (msg) {
        LogError('ERROR Thread $i: $msg');
        tasksRemaining--;
        if (!completer.isCompleted) {
          completer.complete(null);
        }
      });
    }

    return completer.future;
  }

  Future<List<int>> loadFile(String path) {
    LogDebug('REQUEST FILE $path');

    Completer<List<int>> c = new Completer<List<int>>();

    _loadFile(path).then((bytes) {
      LogDebug('LOADED FILE $path');
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
          c.complete(new Uint8List.fromList(bytes.codeUnits));
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
      path = '$scenePath/$path';

      HttpRequest.request(path, method: 'GET',
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
