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
import 'package:dartray/dartray.dart';
import 'package:image/image.dart';

//import 'scenes/01_bowl_of_spheres.dart';
//import 'scenes/02_kon_tiki_directlighting.dart';
//import 'scenes/03_quadrics_directlighting.dart';
//import 'scenes/04_box.dart';
//import 'scenes/05_distant_light.dart';
//import 'scenes/05_distant_light2.dart';
//import 'scenes/06_woman.dart';
//import 'scenes/07_area_light.dart';
//import 'scenes/08_whitted.dart';
//import 'scenes/09_quadrics.dart';

//import 'scenes/bunny.dart';
import 'scenes/cornell_path.dart';
//import 'scenes/area_light.dart';
//import 'scenes/spheres.dart';
//import 'scenes/teapot.dart';
//import 'scenes/room_path.dart';


void main() {
  const int width = 256;
  const int height = 256;

  var c = new Html.CanvasElement();
  Html.document.body.append(c);
  c.width = width;
  c.height = height;

  var imageData = c.context2D.getImageData(0, 0, c.width, c.height);
  var img = new Image(c.width, c.height);


  // Use the custom random number generator so we can verify results against
  // C++ pbrt.
  RNG.UseMathRandom = false;
  //Spectrum.type = Spectrum.SPECTRUM;

  Stopwatch timer = new Stopwatch();
  timer.start();
  new RenderManager().render(SCENE, image: img,
      isolate: 'render_isolate.dart', numThreads: 2,
      log: (int type, String msg) {
        print(msg);
        var div = new Html.Element.html('<div>$msg</div>');
        Html.document.body.nodes.add(div);
      },
      preview: (Image img) {
        var bytes = img.getBytes();
        imageData.data.setRange(0, bytes.length, bytes);
        c.context2D.putImageData(imageData, 0, 0);
      }).then((OutputImage output) {
        timer.stop();
        LogInfo('RENDER FINISHED: ${timer.elapsedMilliseconds / 1000.0} seconds');
        var bytes = img.getBytes();
        imageData.data.setRange(0, bytes.length, bytes);
        c.context2D.putImageData(imageData, 0, 0);
      });
}
