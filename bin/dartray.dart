import 'dart:io';
import 'package:args/args.dart';
import 'package:dartray/dartray.dart';
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

  String scene = new File(args.rest[0]).readAsStringSync();

  Stopwatch timer = new Stopwatch();
  timer.start();
  new RenderManager().render(scene).then((output) {
    timer.stop();
    LogInfo('RENDER FINISHED: ${timer.elapsedMilliseconds / 1000.0} seconds');
    if (output != null) {
      Image image = output.toImage();
      List<int> png = encodeNamedImage(image, out);
      new File(out).writeAsBytesSync(png);
    }
  });
}
