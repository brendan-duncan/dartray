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
import 'dart:io';
import 'package:args/args.dart';
import 'package:dartray/dartray_io.dart';
import 'package:image/image.dart';

void main(List<String> argv) {
  var parser = new ArgParser();
  parser.addOption('output', abbr: 'o', defaultsTo: 'output.png');
  var args = parser.parse(argv);

  if (args.rest.isEmpty) {
    print('Usage: dartray [options] <scene.pbrt>');
    print('options:');
    print(parser.getUsage());
    return;
  }

  String out = args['output'];

  Stopwatch timer = new Stopwatch();
  timer.start();
  new RenderManager().render(args.rest[0]).then((output) {
    timer.stop();
    LogInfo('RENDER FINISHED: ${timer.elapsed}');
    if (output != null) {
      Image image = output.toImage();
      List<int> png = encodeNamedImage(image, out);
      new File(out).writeAsBytesSync(png);
    }
  });
}
