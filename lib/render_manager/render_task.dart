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
 * Manages [RenderIsolate]s controlled by the [RenderManager].
 */
class RenderTask {
  static const int CONNECTING = 1;
  static const int CONNECTED = 2;
  static const int STOPPED = 3;

  int status = CONNECTING;
  ReceivePort receivePort = new ReceivePort();
  SendPort sendPort;
  PreviewCallback previewCallback;
  int taskNum;
  int taskCount;
  Image threadImage;
  List<int> extents = [0, 0, 0, 0];

  RenderTask(this.previewCallback, this.taskNum, this.taskCount);

  Future<OutputImage> render(String scene, Image image, String isolateUri) {
    Completer<OutputImage> completer = new Completer<OutputImage>();

    Sampler.ComputeSubWindow(image.width, image.height, taskNum, taskCount,
                             extents);

    threadImage = new Image.from(image);
    Isolate.spawnUri(Uri.parse(isolateUri), ['_'],
                     receivePort.sendPort).then((iso) {
    });

    receivePort.listen((msg) {
      if (status == CONNECTING) {
        _linkEstablish(msg, image, scene);
      } else if (status == CONNECTED) {
        if (msg is Map && msg.containsKey('cmd')) {
          var cmd = msg['cmd'];

          if (cmd == 'preview' && msg.containsKey('image')) {
            var bytes = msg['image'];
            if (taskCount > 1) {
              threadImage.getBytes().setRange(0, bytes.length, bytes);
              copyInto(image, threadImage, dstX: extents[0], dstY: extents[2],
                       srcX: extents[0], srcY: extents[2],
                       srcW: extents[1], srcH: extents[3]);
            } else {
              image.getBytes().setRange(0, bytes.length, bytes);
            }
            if (previewCallback != null) {
              previewCallback(image);
            }
            return;
          } else if (cmd == 'error') {
            LogInfo('ERROR: ${msg['msg']}');
            completer.completeError(msg['msg']);
            return;
          } else if (cmd == 'final' && msg.containsKey('output')) {
            Float32List rgb = msg['output'];
            OutputImage output = new OutputImage(extents[0], extents[2],
                                                 extents[1], extents[3],
                                                 rgb);
            /*if (taskCount > 1) {
              threadImage.getBytes().setRange(0, bytes.length, bytes);
              copyInto(image, threadImage, dstX: extents[0], dstY: extents[2],
                       srcX: extents[0], srcY: extents[2],
                       srcW: extents[1], srcH: extents[3]);
            } else {
              image.getBytes().setRange(0, bytes.length, bytes);
            }

            if (previewCallback != null) {
              previewCallback(image);
            }*/

            completer.complete(output);
            return;
          }
        }

        LogInfo(msg);
      }
    });

    return completer.future;
  }

  void _linkEstablish(msg, [Image image, String scene]) {
    if (msg is SendPort) {
      sendPort = msg;
      sendPort.send('ping');
    } else if (msg == 'pong') {
      status = CONNECTED;
      _startIsolateRender(image, scene);
    }
  }

  void _startIsolateRender(Image image, String scene) {
    sendPort.send({'cmd': 'render', 'width': image.width,
                   'height': image.height, 'scene': scene,
                   'taskNum': taskNum, 'taskCount': taskCount,
                   'preview': previewCallback != null ? true : false});
  }
}
