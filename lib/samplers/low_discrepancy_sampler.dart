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

class LowDiscrepancySampler extends Sampler {
  LowDiscrepancySampler(int x, int y, int width, int height,
                        double sopen, double sclose, this.pixels,
                        int nsamp) :
    super(x, y, width, height, sopen, sclose, RoundUpPow2(nsamp)) {
    if (pixels == null) {
      LogSevere('A PixelSampler is required by LowDiscrepencySampler');
    }

    pixels.setup(x, y, width, height);
    pixelIndex = 0;

    if (!IsPowerOf2(nsamp)) {
      nPixelSamples = RoundUpPow2(nsamp);
      LogWarning('Pixel samples being rounded up to power of 2: '
                 '$nsamp => $nPixelSamples');
    } else {
      nPixelSamples = nsamp;
    }
    sampleBuf = null;

    pass = 0;
    if (RenderOverrides.SamplingMode() == Sampler.TWO_PASS_SAMPLING ||
        RenderOverrides.SamplingMode() == Sampler.ITERATIVE_SAMPLING) {
      randomSampler = new RandomSampler(x, y, width, height,
                                        sopen, sclose, pixels, 1);
    }
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

  static LowDiscrepancySampler Create(ParamSet params, int x, int y, int width,
                                      int height, Camera camera,
                                      PixelSampler pixels) {
    // Initialize common sampler parameters
    int nsamp = params.findOneInt('pixelsamples', 4);

    return new LowDiscrepancySampler(x, y, width, height, camera.shutterOpen,
                                     camera.shutterClose, pixels, nsamp);
  }

  PixelSampler pixels;
  Int32List pixel = new Int32List(2);
  int pixelIndex;
  int nPixelSamples;
  Float32List sampleBuf;
  Sampler randomSampler;
  int pass;
}
