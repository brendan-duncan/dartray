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
part of render_manager;

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

  void start(SendPort port) {
    Log = _log;
    status = CONNECTING;
    receivePort = new ReceivePort();
    sendPort = port;

    sendPort.send(receivePort.sendPort);

    receivePort.listen((msg) {
      if (status == CONNECTING) {
        _linkEstablish(msg);
      } else if (status == CONNECTED) {
        _run(msg);
      }
    });
  }

  void _linkEstablish(msg) {
    if (msg == 'ping') {
      sendPort.send('pong');
      status = CONNECTED;
    }
  }

  void _run([msg = null]) {
    if (msg is Map) {
      if (msg.containsKey('cmd')) {
        var cmd = msg['cmd'];
        if (cmd == 'render') {
          taskNum = msg.containsKey('taskNum') ? msg['taskNum'] : 0;
          taskCount = msg.containsKey('taskCount') ? msg['taskCount'] : 1;
          int w = msg.containsKey('width') ? msg['width'] : 256;
          int h = msg.containsKey('height') ? msg['height'] : 256;
          String scene = msg.containsKey('scene') ? msg['scene'] : '';
          bool doPreview = msg.containsKey('preview') ? msg['preview'] : false;

          _render(w, h, scene, taskNum, taskCount, doPreview);
        }
      }
    } else if (msg == 'quit') {
      status = STOPPED;
      receivePort.close();
    }
  }

  void _log(int type, String msg) {
    String timestamp = new DateTime.now().toString().substring(11);
    sendPort.send('${LOG_TYPES[type]} [THREAD $taskNum/$taskCount]: '
                  '$timestamp : $msg');
  }

  bool _render(int w, int h, String scene, int taskNum, int taskCount,
               bool doPreview) {
    _log(LOG_INFO, 'RENDER ${w}x${h}');
    Stopwatch timer = new Stopwatch()..start();

    var img = new Image(w, h);
    Pbrt pbrt = new Pbrt();

    if (doPreview) {
      pbrt.setPreviewCallback((Image img) {
        sendPort.send({'cmd': 'preview', 'image': img.getBytes()});
      });
    }

    OutputImage output;
    try {
      pbrt.setTask(taskNum, taskCount);
      output = pbrt.renderScene(scene, img);
    } catch (e) {
      sendPort.send({'cmd': 'error', 'msg': e.toString()});
      return false;
    }

    _log(LOG_INFO, 'FINISHED: ${timer.elapsed}');
    LogInfo('[$taskNum] Stats....\n${Stats.getString()}');

    sendPort.send({'cmd': 'final', 'output': output.rgb});

    return true;
  }
}
