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
part of samplers;

class RandomSampler extends Sampler {
  RandomSampler(int x, int y, int width, int height, double sopen,
                double sclose, this.pixels, int ns) :
    super(x, y, width, height, sopen, sclose, ns) {
    if (pixels == null) {
      LogSevere('A PixelSampler is required by RandomSampler');
    }
    pixels.setup(x, y, width, height);
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
    }

    pixels.getPixel(pixelIndex++, pixel);

    int mode = RenderOverrides.SamplingMode();
    int numSamples = mode == Sampler.ITERATIVE_SAMPLING ? 1 :
                     mode == Sampler.FULL_SAMPLING ? samplesPerPixel :
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

  static RandomSampler Create(ParamSet params, int x, int y, int width,
                              int height, Camera camera, PixelSampler pixels) {
    int ns = params.findOneInt('pixelsamples', 10);

    return new RandomSampler(x, y, width, height,
                             camera.shutterOpen, camera.shutterClose,
                             pixels, ns);
  }

  PixelSampler pixels;
  Int32List pixel = new Int32List(2);
  int pixelIndex;
  int sampleCount;
}
