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

class LowDiscrepancySampler extends Sampler {
  LowDiscrepancySampler(int xstart, int xend, int ystart, int yend,
                        double sopen, double sclose, this.pixels,
                        int nsamp, int samplingMode) :
    super(xstart, xend, ystart, yend, sopen, sclose,
          RoundUpPow2(nsamp), samplingMode) {
    if (pixels == null) {
      LogSevere('A PixelSampler is required by LowDiscrepencySampler');
    }
    pixels.setup(xstart, xend, ystart, yend);
    pixelIndex = 0;
    if (!IsPowerOf2(nsamp)) {
      LogWarning('Pixel samples being rounded up to power of 2');
      nPixelSamples = RoundUpPow2(nsamp);
    } else {
      nPixelSamples = nsamp;
    }
    sampleBuf = null;

    if (samplingMode == Sampler.TWO_PASS_SAMPLING ||
        samplingMode == Sampler.ITERATIVE_SAMPLING) {
      randomSampler = new RandomSampler(xstart, xend, ystart, yend,
                                        sopen, sclose, pixels, 1,
                                        Sampler.ITERATIVE_SAMPLING);
    }

    pass = 0;
  }

  Sampler getSubSampler(int num, int count) {
    List extents = [0, 0, 0, 0];
    computeSubWindow(num, count, extents);
    if (extents[0] == extents[1] || extents[2] == extents[3]) {
      return null;
    }

    return new LowDiscrepancySampler(extents[0], extents[1],
                                     extents[2], extents[3],
                                     shutterOpen, shutterClose,
                                     pixels, nPixelSamples, samplingMode);
  }

  int roundSize(int size) {
    return RoundUpPow2(size);
  }

  int getMoreSamples(List<Sample> samples, RNG rng) {
    if (pass == 0 && randomSampler != null) {
      int count = randomSampler.getMoreSamples(samples, rng);
      if (count != 0) {
        return count;
      }
      pass++;
    }

    if (pixelIndex >= pixels.numPixels()) {
      return 0;
    }

    if (sampleBuf == null) {
       sampleBuf = new Float32List(LDPixelSampleFloatsNeeded(samples[0],
                                                            nPixelSamples));
    }

    pixels.getPixel(pixelIndex++, pixel);

    LDPixelSample(pixel[0], pixel[1], shutterOpen, shutterClose,
                 nPixelSamples, samples, sampleBuf, rng);

    return nPixelSamples;
  }

  int maximumSampleCount() {
    return nPixelSamples;
  }

  static LowDiscrepancySampler Create(ParamSet params, Film film,
                                      Camera camera, PixelSampler pixels) {
    // Initialize common sampler parameters
    List<int> extents = [0, 0, 0, 0];
    film.getSampleExtent(extents);
    int nsamp = params.findOneInt('pixelsamples', 4);

    String mode = params.findOneString('mode', 'full');
    int samplingMode = (mode == 'full') ? Sampler.FULL_SAMPLING :
                       (mode == 'twopass') ? Sampler.TWO_PASS_SAMPLING :
                       (mode == 'iterative') ? Sampler.ITERATIVE_SAMPLING :
                       -1;
    if (samplingMode == -1) {
      LogWarning('Invalid sampling mode: $mode. Using \'full\'.');
      samplingMode = Sampler.FULL_SAMPLING;
    }

    return new LowDiscrepancySampler(extents[0], extents[1], extents[2],
                                     extents[3], camera.shutterOpen,
                                     camera.shutterClose, pixels,
                                     nsamp, samplingMode);
  }

  PixelSampler pixels;
  Int32List pixel = new Int32List(2);
  int pixelIndex;
  int nPixelSamples;
  Float32List sampleBuf;
  Sampler randomSampler;
  int pass;
}
