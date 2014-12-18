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
import 'dart:html';
import 'package:crypto/crypto.dart';
import 'package:dartray/dartray_web.dart';
import 'package:image/image.dart';

const List<String> SCENES = const [
  'anim-bluespheres',
  'bump-sphere',
  'bunny',
  'cornell-mlt',
  'cornell-path',
  'miscquads',
  'smoke-2',
  'spotfog',
  'teapot-area-light',
  'teapot-metal'
];

void main() {
  String scene = 'cornell-path';

  String queryString = window.location.search.replaceFirst('?', '');
  if (queryString.isNotEmpty) {
    scene = queryString;
  }

  if (!SCENES.contains(scene)) {
    LogError('Unknown Scene');
    return;
  }

  var sceneMenu = querySelector('#sceneMenu');
  for (var s in SCENES) {
    var item = new OptionElement();
    item.value = s;
    item.text = s;
    if (s == scene) {
      item.selected = true;
    }
    sceneMenu.append(item);
  }
  sceneMenu.onChange.listen((e) {
    String s = SCENES[sceneMenu.selectedIndex];
    var url = '${window.location.origin}/${window.location.pathname}?${s}';
    window.location.replace(url);
  });

  scene = 'scenes/$scene.pbrt';

  RenderManager renderManager = new RenderManager();

  var log = querySelector('#log');
  var canvas = querySelector('#renderCanvas');
  var context = canvas.context2D;
  var imageData = context.getImageData(0, 0, canvas.width, canvas.height);
  var canvasContainer = querySelector('#canvasContainer');

  RenderOverrides overrides = new RenderOverrides();
  overrides.samplingMode = Sampler.TWO_PASS_SAMPLING;

  Stopwatch timer = new Stopwatch();
  timer.start();
  renderManager.render(scene,
      isolate: 'dartray_isolate.dart',
      //numThreads: 2,
      overrides: overrides,
      log: (int type, String msg) {
        print('$msg');
        if (type != LOG_DEBUG) {
          log.text += '$msg\n';
          // TODO this seems to be missing from preElement now...
          //log.scrollByLines(1);
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
          // When the render has finished, replace the canvas with an IMG
          // element so it can be saved from the browser.
          var img = output.toImage(gamma: 2.2);
          var png = encodePng(img);

          var imgEl = new ImageElement();
          imgEl.id = 'renderImage';
          canvasContainer.append(imgEl);

          canvas.hidden = true;

          var base64 = CryptoUtils.bytesToBase64(png);
          imgEl.src = 'data:image/png;base64,${base64}';
        }
      });
}
