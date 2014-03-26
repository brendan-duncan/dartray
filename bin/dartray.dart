import 'dart:io';
import 'package:args/args.dart';
import 'package:dartray/dartray.dart';
import 'package:image/image.dart';

void main(List<String> argv) {
  var parser = new ArgParser();
  var args = parser.parse(argv);

  if (args.rest.isEmpty) {
    print('Usage: dartray <scene.pbrt>');
  }

  String scene = new File(args.rest[0]).readAsStringSync();

  Stopwatch timer = new Stopwatch();
  timer.start();
  new RenderManager().render(scene, null, write: writeImage).then((e) {
    timer.stop();
    LogInfo('RENDER FINISHED: ${timer.elapsedMilliseconds / 1000.0} seconds');
  });
}

void writeImage(String name, List<double> pixels, int xRes,
                int yRes, int totalXRes, int totalYRes,
                int xOffset, int yOffset) {
  print('WRITE $name');
}
