part of core;

/**
 * Storage for the output of a renderer.
 *
 * A render thread may work on a portion of the overal image. Only the output
 * of the threads portion of the overal image is stored. The [xOffset],
 * [yOffset], [width], and [height] parameters specify the region of the
 * overal image it occupies.
 *
 * Pixel colors are stored in high dynamic range floating-point values.
 * Tone-mapping will need to be done to translate the pixels into a low dynamic
 * range format.
 */
class OutputImage {
  int width;
  int height;
  int xOffset;
  int yOffset;
  final Float32List rgb;

  OutputImage(this.xOffset, this.yOffset, int width, height, [Float32List rgb]) :
    this.width = width,
    this.height = height,
    this.rgb = rgb != null ? rgb : new Float32List(width * height * 3);

  Img.Image toImage({double gamma: 2.2}) {
    Img.Image img = new Img.Image(width, height);
    Uint8List pixels = img.getBytes();
    for (int i = 0, oi = 0, len = rgb.length; i < len; i += 3, oi += 4) {
      pixels[oi] = (rgb[i] * 255.0).floor().clamp(0, 255);
      pixels[oi + 1] = (rgb[i + 1] * 255.0).floor().clamp(0, 255);
      pixels[oi + 2] = (rgb[i + 2] * 255.0).floor().clamp(0, 255);
      pixels[oi + 3] = 255;
    }
    return Img.adjustColor(img, gamma: 1.0/gamma);
  }
}
