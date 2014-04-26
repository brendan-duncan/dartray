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

/**
 * Starts an Isolate service to connect to the primary [RenderManager] from
 * a Dart Isolate.
 *
 * A RenderIsolate can render a portion of the image, such that multiple
 * RenderIsolates can run in parallel, providing multi-threading for the
 * renderer.
 *
 * Some complications of Isolates:
 * 1) There is no shared memory between isolates or the main thread. That
 * means every isolate needs its own copy of the scene; of textures, geometry
 * and other resources; acceleration structures and other caches. For large
 * scenes this could be problematic.
 *
 * 2) All data shared between isolates happens through messages sent over a
 * port. Only sting or byte buffers can be sent (maps and lists can be sent
 * too, where they are serialized prior), and no shared memory so these buffers
 * are copied. Therefore, all communication between isolates have significant
 * overhead.
 */
class RenderIsolate {
  final CONNECTING = 1;
  final CONNECTED = 2;
  final STOPPED = 3;

  int status;
  int taskNum = 0;
  int taskCount = 1;
  ReceivePort receivePort;
  SendPort sendPort;
  RenderManagerInterface manager;
  Math.Random rng = new Math.Random();

  RenderIsolate(this.manager);

  void start(SendPort port) {
    Log = _log;
    status = CONNECTING;
    receivePort = new ReceivePort();
    sendPort = port;

    sendPort.send(receivePort.sendPort);
    status = CONNECTED;

    receivePort.listen((msg) {
      _message(msg);
    });
  }

  /**
   * Sends a request to the parent isolate, sending a unique id.
   * The returned [Future] will be completed when the parent isolate
   * sends a message back with that ID.
   */
  Future requestResponse(msg) {
    int id = rng.nextInt(0xffffffff);
    Completer c = new Completer();
    Map cmd = {'cmd': 'request', 'id': id, 'msg': msg};
    sendPort.send(cmd);
    requests[id] = c;
    return c.future;
  }

  /**
   * Handle a message recieved from the parent isolate.
   */
  void _message([msg = null]) {
    if (msg is Map) {
      if (msg.containsKey('cmd')) {
        var cmd = msg['cmd'];
        if (cmd == 'request') {
          int id = msg['id'];
          if (requests.containsKey(id)) {
            Completer c = requests[id];
            requests.remove(id);
            c.complete(msg['data']);
          }
        } else if (cmd == 'render') {
          taskNum = msg.containsKey('taskNum') ? msg['taskNum'] : 0;
          taskCount = msg.containsKey('taskCount') ? msg['taskCount'] : 1;
          String scene = msg.containsKey('scene') ? msg['scene'] : '';
          bool doPreview = msg.containsKey('preview') ? msg['preview'] : false;

          RenderOverrides overrides;
          try {
            overrides = msg.containsKey('overrides') ?
                        new RenderOverrides.fromJson(msg['overrides']) :
                        null;
          } catch (e) {
            overrides = null;
          }

          _render(scene, taskNum, taskCount, doPreview, overrides);
        }
      }
    } else if (msg == 'pause') {
      LogInfo('Isolate $taskNum/$taskCount PAUSE');
    } else if (msg == 'resume') {
      LogInfo('Isolate $taskNum/$taskCount RESUME');
    } else if (msg == 'stop') {
      LogInfo('Isolate $taskNum/$taskCount STOP');
      status = STOPPED;
      receivePort.close();
    }
  }

  void _log(int type, String msg) {
    String timestamp = new DateTime.now().toString().substring(11);
    sendPort.send('${LOG_TYPES[type]} [THREAD ${taskNum + 1}/$taskCount]: '
                  '$timestamp : $msg');
  }

  void _render(String scene, int taskNum, int taskCount,
               bool doPreview, RenderOverrides overrides) {
    LogInfo('RENDER THREAD STARTED ${taskNum + 1} / $taskCount');

    Stopwatch timer = new Stopwatch()..start();

    DartRay dartray = new DartRay(manager);

    List<int> extents = [0, 0, 0, 0];

    if (doPreview) {
      dartray.setPreviewCallback((Image img) {
        GetSubWindow(img.width, img.height, taskNum, taskCount, extents);

        sendPort.send({'cmd': 'preview',
                       'res': [img.width, img.height],
                       'extents': extents,
                       'image': img.getBytes()});
      });
    }

    try {
      OutputImage output;
      dartray.setTask(taskNum, taskCount);
      dartray.renderScene(scene, overrides: overrides).then((output) {
        _log(LOG_INFO, 'FINISHED: ${timer.elapsed}');
        LogInfo('[$taskNum] STATS:\n${Stats.getString()}');

        GetSubWindow(output.width, output.height, taskNum, taskCount, extents);

        sendPort.send({'cmd': 'final',
                       'output': output.rgb,
                       'res': [output.width, output.height],
                       'extents': extents});
      });
    } catch (e) {
      LogError('ERROR: ${e}');
      sendPort.send({'cmd': 'error', 'msg': e.toString()});
    }
  }

  Map<int, Completer> requests = {};
}
