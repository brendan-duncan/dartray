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
part of samplers;

class RandomSampler extends Sampler {
  RandomSampler(int xstart, int xend, int ystart,
      int yend, double sopen, double sclose, this.pixels,
      int ns, int samplingMode) :
    super(xstart, xend, ystart, yend, sopen, sclose, ns, samplingMode) {
    if (pixels == null) {
      LogSevere('A PixelSampler is required by RandomSampler');
    }
    pixels.setup(xstart, xend, ystart, yend);
    pixelIndex = 0;
    sampleCount = 0;
  }

  int maximumSampleCount() {
    return samplesPerPixel;
  }

  int getMoreSamples(List<Sample> sample, RNG rng) {
    PixelSampler pixels = this.pixels;

    if (pixelIndex >= pixels.numPixels()) {
      sampleCount++;
      if (sampleCount >= samplesPerPixel) {
        return 0;
      }
      pixelIndex = 0;
      LogInfo('Random Sampler PASS $sampleCount');
    }

    pixels.getPixel(pixelIndex++, pixel);

    int numSamples = samplingMode == Sampler.ITERATIVE_SAMPLING ? 1 :
                     samplingMode == Sampler.FULL_SAMPLING ? samplesPerPixel :
                     (sampleCount == 0) ? 1 : samplesPerPixel - 1;

    for (int si = 0; si < numSamples; ++si) {
      // Return next sample point
      sample[si].imageX = rng.randomFloat() + pixel[0];
      sample[si].imageY = rng.randomFloat() + pixel[1];
      sample[si].lensU = rng.randomFloat();
      sample[si].lensV = rng.randomFloat();
      sample[si].time = Lerp(rng.randomFloat(), shutterOpen, shutterClose);

      // Generate stratified samples for integrators
      for (int i = 0; i < sample[si].n1D.length; ++i) {
        for (int j = 0; j < sample[si].n1D[i]; ++j) {
          sample[si].oneD[i][j] = rng.randomFloat();
        }
      }

      for (int i = 0; i < sample[si].n2D.length; ++i) {
        for (int j = 0; j < 2 * sample[si].n2D[i]; ++j) {
          sample[si].twoD[i][j] = rng.randomFloat();
        }
      }
    }

    return numSamples;
  }

  int roundSize(int sz) {
    return sz;
  }

  Sampler getSubSampler(int num, int count) {
    List<int> extents = [0, 0, 0, 0];
    computeSubWindow(num, count, extents);
    if (extents[0] == extents[1] || extents[2] == extents[3]) {
      return null;
    }

    return new RandomSampler(extents[0], extents[1], extents[2], extents[3],
                             shutterOpen, shutterClose, pixels,
                             samplesPerPixel, samplingMode);
  }

  static RandomSampler Create(ParamSet params, Film film, Camera camera,
                              PixelSampler pixels) {
    int ns = params.findOneInt('pixelsamples', 32);

    List<int> extents = [0, 0, 0, 0];
    film.getSampleExtent(extents);

    String mode = params.findOneString('mode', 'full');
    int samplingMode = (mode == 'full') ? Sampler.FULL_SAMPLING :
                       (mode == 'twopass') ? Sampler.TWO_PASS_SAMPLING :
                       (mode == 'iterative') ? Sampler.ITERATIVE_SAMPLING :
                       -1;
    if (samplingMode == -1) {
      LogWarning('Invalid sampling mode: $mode. Using \'full\'.');
      samplingMode = Sampler.FULL_SAMPLING;
    }
    LogInfo(mode);

    return new RandomSampler(extents[0], extents[1], extents[2], extents[3],
                             camera.shutterOpen, camera.shutterClose,
                             pixels, ns, samplingMode);
  }

  PixelSampler pixels;
  Int32List pixel = new Int32List(2);
  int pixelIndex;
  int sampleCount;
}
