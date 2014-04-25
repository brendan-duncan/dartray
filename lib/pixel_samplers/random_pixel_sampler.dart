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
part of pixel_samplers;

/**
 * Sample the image pixels in a random order.
 */
class RandomPixelSampler extends PixelSampler {
  void setup(int x, int y, int width, int height) {
    super.setup(x, y, width, height);

    _samples = new Int32List(width * height * 2);

    for (int y = top, si = 0; y <= bottom; ++y) {
      for (int x = left; x <= right; ++x) {
        _samples[si++] = x;
        _samples[si++] = y;
      }
    }

    RNG rng = new RNG();

    _numSamples = _samples.length ~/ 2;
    // Shuffle the samples
    for (int i = 0, r = 0; i < _numSamples; ++i, r += 2) {
      int l = (rng.randomUint() % _numSamples) * 2;

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

  static RandomPixelSampler Create(ParamSet params) {
    return new RandomPixelSampler();
  }

  int _numSamples;
  Int32List _samples;
}
