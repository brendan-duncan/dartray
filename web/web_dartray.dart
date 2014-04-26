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
import 'dart:html' as Html;
import 'package:crypto/crypto.dart';
import 'package:dartray/dartray_web.dart';
import 'package:image/image.dart';

//String scene = 'scenes/bowl_of_spheres.pbrt';
//String scene = 'scenes/kon_tiki_directlighting.pbrt';
//String scene = 'scenes/quadrics_directlighting.pbrt';
//String scene = 'scenes/box.pbrt';
//String scene = 'scenes/distant_light.pbrt';
//String scene = 'scenes/whitted.pbrt';
//String scene = 'scenes/quadrics.pbrt';
//String scene = 'scenes/area_light.pbrt';
//String scene = 'scenes/bunny.pbrt';
String scene = 'scenes/cornell_path.pbrt';
//String scene = 'scenes/room-path.pbrt';
//String scene = 'scenes/spheres.pbrt';
//String scene = 'scenes/teapot-area-light.pbrt';
//String scene = 'scenes/nurbs.pbrt';

void main() {
  RenderManager renderManager = new RenderManager();

  var log = Html.querySelector('#log');

  var canvas = Html.querySelector('#renderCanvas');
  var context = canvas.context2D;
  var imageData = context.getImageData(0, 0, canvas.width, canvas.height);

  var canvasContainer = Html.querySelector('#canvasContainer');

  RenderOverrides overrides = new RenderOverrides();
  //overrides.setSampler('random', {'integer pixelsamples': [1]});
  overrides.resolutionScale = 0.5;
  overrides.samplingMode = Sampler.TWO_PASS_SAMPLING;

  Stopwatch timer = new Stopwatch();
  timer.start();
  renderManager.render(scene,
      isolate: 'web_isolate.dart',
      numThreads: 2,
      overrides: overrides,
      log: (int type, String msg) {
        print('$msg');
        if (type != LOG_DEBUG) {
          log.text += '$msg\n';
          log.scrollByLines(1);
        }
      },
      preview: (Image img) {
        if (img.width != canvas.width || img.height != canvas.height) {
          canvas.width = img.width;
          canvas.height = img.height;
          imageData = context.getImageData(0, 0, canvas.width, canvas.height);
        }
        var bytes = img.getBytes();
        imageData.data.setRange(0, bytes.length, bytes);
        context.putImageData(imageData, 0, 0);
      }).then((OutputImage output) {
        timer.stop();
        LogInfo('FINISHED Render: ${timer.elapsed}');

        String s = Stats.getString();
        if (s.isNotEmpty) {
          LogInfo('STATS....\n${Stats.getString()}');
        }

        if (output != null) {
          Image img = output.toImage(gamma: 2.2);
          if (img.width != canvas.width || img.height != canvas.height) {
            canvas.width = img.width;
            canvas.height = img.height;
            imageData = context.getImageData(0, 0, canvas.width, canvas.height);
          }

          // When the render has finished, replace the canvas with an IMG
          // element so it can be saved from the browser.
          var png = encodePng(img);

          var imgEl = new Html.ImageElement();
          imgEl.id = 'renderImage';
          canvasContainer.append(imgEl);

          canvas.hidden = true;

          var png64 = CryptoUtils.bytesToBase64(png);
          imgEl.src = 'data:image/png;base64,${png64}';
        }
      });
}
