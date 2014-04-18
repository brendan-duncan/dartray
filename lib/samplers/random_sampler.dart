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
  static RandomSampler Create(ParamSet params, Film film, Camera camera,
                              PixelSampler pixels) {
    int ns = params.findOneInt('pixelsamples', 10);
    bool continuous = params.findOneBool('continuous', false);
    List<int> extents = [0, 0, 0, 0];
    film.getSampleExtent(extents);
    return new RandomSampler(extents[0], extents[1], extents[2], extents[3],
                             ns, continuous,
                             camera.shutterOpen, camera.shutterClose,
                             pixels);
  }

  RandomSampler(int xstart, int xend, int ystart,
      int yend, int ns, this.continuous, double sopen, double sclose,
      this.pixels) :
    super(xstart, xend, ystart, yend, ns, sopen, sclose) {
    if (pixels == null) {
      LogSevere('Pixel sampler is required by LowDiscrepencySampler');
    }
    pixels.setup(xstart, xend, ystart, yend);
    pixelIndex = 0;
    sampleCount = 0;
    // Get storage for a pixel's worth of stratified samples
    imageSamples = new Float32List(2);
    lensSamples = new Float32List(2);
    timeSamples = new Float32List(1);
  }

  int maximumSampleCount() {
    return 1;
  }

  int getMoreSamples(List<Sample> sample, RNG rng) {
    if (pixelIndex >= pixels.numPixels()) {
      sampleCount++;
      if (!continuous) {
        if (sampleCount >= samplesPerPixel) {
          return 0;
        }
      }
      pixelIndex = 0;
    }

    pixels.getPixel(pixelIndex++, pixel);

    timeSamples[0] = rng.randomFloat();
    imageSamples[0] = rng.randomFloat();
    lensSamples[0] = rng.randomFloat();

    imageSamples[1] = rng.randomFloat();
    lensSamples[1] = rng.randomFloat();

    // Shift image samples to pixel coordinates
    imageSamples[0] += pixel[0];
    imageSamples[1] += pixel[1];

    // Return next sample point
    sample[0].imageX = imageSamples[0];
    sample[0].imageY = imageSamples[1];
    sample[0].lensU = lensSamples[0];
    sample[0].lensV = lensSamples[1];
    sample[0].time = Lerp(timeSamples[0], shutterOpen, shutterClose);

    // Generate stratified samples for integrators
    for (int i = 0; i < sample[0].n1D.length; ++i) {
      for (int j = 0; j < sample[0].n1D[i]; ++j) {
        sample[0].oneD[i][j] = rng.randomFloat();
      }
    }

    for (int i = 0; i < sample[0].n2D.length; ++i) {
      for (int j = 0; j < 2 * sample[0].n2D[i]; ++j) {
        sample[0].twoD[i][j] = rng.randomFloat();
      }
    }

    return 1;
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
                             samplesPerPixel, continuous, shutterOpen,
                             shutterClose, pixels);
  }

  bool continuous;
  PixelSampler pixels;
  Int32List pixel = new Int32List(2);
  int pixelIndex;
  Float32List imageSamples;
  Float32List lensSamples;
  Float32List timeSamples;
  int sampleCount;
}
