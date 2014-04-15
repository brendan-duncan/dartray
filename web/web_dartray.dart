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
import 'dart:html' as Html;
import 'package:dartray/dartray_web.dart';
import 'package:image/image.dart';

//String scene = 'scenes/01_bowl_of_spheres.pbrt';
//String scene = 'scenes/02_kin_tiki_directlighting.pbrt';
//String scene = 'scenes/03_quadrics_directlighting.pbrt';
//String scene = 'scenes/04_box.pbrt';
//String scene = 'scenes/05_distant_light.pbrt';
//String scene = 'scenes/07_area_light.pbrt';
//String scene = 'scenes/08_whitted.pbrt';
//String scene = 'scenes/09_quadrics.pbrt';
//String scene = 'scenes/area_light.pbrt';
//String scene = 'scenes/bunny.pbrt';
String scene = 'scenes/cornell_path.pbrt';
//String scene = 'scenes/room-path.pbrt';
//String scene = 'scenes/spheres.pbrt';
//String scene = 'scenes/teapot-area-light.pbrt';

void main() {
  const int width = 256;
  const int height = 256;

  var c = new Html.CanvasElement();
  Html.document.body.append(c);
  c.width = width;
  c.height = height;

  c.context2D.fill();
  var imageData = c.context2D.getImageData(0, 0, c.width, c.height);
  var img = new Image(c.width, c.height);

  Stopwatch timer = new Stopwatch();
  timer.start();
  new RenderManager().render(scene,
      image: img,
      isolate: 'web_isolate.dart',
      log: (int type, String msg) {
        print('$msg');
        var div = new Html.Element.html('<pre>$msg</pre>');
        Html.document.body.nodes.add(div);
      },
      preview: (Image img) {
        var bytes = img.getBytes();
        imageData.data.setRange(0, bytes.length, bytes);
        c.context2D.putImageData(imageData, 0, 0);
      }).then((OutputImage output) {
        timer.stop();
        LogInfo('RENDER FINISHED: ${timer.elapsed}');
        String s = Stats.getString();
        if (s.isNotEmpty) {
          LogInfo('STATS....\n${Stats.getString()}');
        }
        var bytes = img.getBytes();
        imageData.data.setRange(0, bytes.length, bytes);
        c.context2D.putImageData(imageData, 0, 0);
      });
}
