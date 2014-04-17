part of pixel_samplers;

/**
 * Sample the image pixels in a random order.
 */
class RandomPixelSampler extends PixelSampler {
  static RandomPixelSampler Create(ParamSet params, Film film) {
    return new RandomPixelSampler();
  }

  void setup(int xPixelStart, int xPixelEnd, int yPixelStart,
             int yPixelEnd) {
    super.setup(xPixelStart, xPixelEnd, yPixelStart, yPixelEnd);

    _samples = new Int32List((xPixelEnd - xPixelStart) *
                             (yPixelEnd - yPixelStart) * 2);

    for (int y = yPixelStart, si = 0; y < yPixelEnd; ++y) {
      for (int x = xPixelStart; x < xPixelEnd; ++x) {
        _samples[si++] = x;
        _samples[si++] = y;
      }
    }

    RNG rng = new RNG();

    _numSamples = _samples.length ~/ 2;
    // Shuffle the samples
    for (int i = 0, r = 0; i < _numSamples; ++i, r += 2) {
      int l = (rng.randomUInt() % _numSamples) * 2;

      // swap x coordinate
      int t = _samples[r];
      _samples[r] = _samples[l];
      _samples[l] = t;

      // swap y coordinate
      t = _samples[r + 1];
      _samples[r + 1] = _samples[l + 1];
      _samples[l + 1] = t;
    }
  }

  int numPixels() {
    return _numSamples;
  }

  void getPixel(int index, List<int> pixel) {
    index *= 2;
    if (index >= _samples.length - 1) {
      return;
    }
    pixel[0] = _samples[index];
    pixel[1] = _samples[index + 1];
  }

  int _numSamples;
  Int32List _samples;
}
