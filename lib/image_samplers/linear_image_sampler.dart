part of image_samplers;

/**
 * Sample the image pixels in a linear order.
 */
class LinearImageSampler extends ImageSampler {
  LinearImageSampler(int xPixelStart, int xPixelEnd, int yPixelStart,
                     int yPixelEnd) :
    super(xPixelStart, xPixelEnd, yPixelStart, yPixelEnd),
    _samples = new Int32List((xPixelEnd - xPixelStart) *
                             (yPixelEnd - yPixelStart) * 2) {
    for (int y = yPixelStart, si = 0; y < yPixelEnd; ++y) {
      for (int x = xPixelStart; x < xPixelEnd; ++x) {
        _samples[si++] = x;
        _samples[si++] = y;
      }
    }

    _numSamples = _samples.length ~/ 2;
  }

  int numPixels() => _numSamples;

  void getPixel(int index, List<int> pixel) {
    index *= 2;
    if (index >= _samples.length - 1) {
      return;
    }
    pixel[0] = _samples[index];
    pixel[1] = _samples[index + 1];
  }

  int _numSamples;
  final Int32List _samples;
}
