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
library render_manager;

import 'dart:async';
import 'dart:isolate';
import 'dart:typed_data';
import 'package:image/image.dart';

import '../core/core.dart';
import '../pbrt/pbrt.dart';

part 'render_isolate.dart';
part 'render_task.dart';

/**
 * Manages the rendering process, either rendering locally or submitting the
 * job to one or more isolates (web workers).
 */
class RenderManager {
  Pbrt pbrt = new Pbrt();

  RenderManager();

  /**
   * [isolateUri] is the path to the script that initiates the
   * RenderIsolate.  This script usually would do nothing more than:
   * import 'dart:isolate';
   * import 'package:dartray/dartray.dart';
   * void main(List<String> args, SendPort port) {
   *   new RenderIsolate().start(port);
   * }
   */
  Future<OutputImage> render(String scene, {Image image, String isolate,
              LogCallback log, PreviewCallback preview,
              int numThreads: 1}) {
    if (log != null) {
      Log = log;
    }

    Completer<OutputImage> completer = new Completer<OutputImage>();

    if (isolate == null) {
      LogInfo('STARTING RENDER');
      OutputImage output = pbrt.renderScene(scene, image);

      if (preview != null) {
        preview(image);
      }

      LogInfo('....rayIntersection: ${Stats.rayIntersection}');
      LogInfo('....rayIntersectionP: ${Stats.rayIntersectionP}');
      completer.complete(output);

      return completer.future;
    }

    int tasksRemaining = numThreads;
    List<RenderTask> jobs = new List<RenderTask>(numThreads);
    for (int i = 0; i < numThreads; ++i) {
      jobs[i] = new RenderTask(preview, i, numThreads);
      jobs[i].render(scene, image, isolate).then((task) {
        tasksRemaining--;
        if (tasksRemaining == 0) {
          completer.complete();
        }
      }, onError: (msg) {
        LogError('Error Thread $i: $msg');
        tasksRemaining--;
        completer.complete(null);
      });
    }

    return completer.future;
  }
}

