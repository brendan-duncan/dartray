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

class HaltonSampler extends Sampler {
  HaltonSampler(int x, int y, int width, int height, double sopen,
                double sclose, int ps) :
    super(x, y, width, height, sopen, sclose, ps) {
    int delta = Math.max(width, height);
    wantedSamples = samplesPerPixel * delta * delta;
    currentSample = 0;

    pass = 0;
    if (RenderOverrides.SamplingMode() == Sampler.TWO_PASS_SAMPLING ||
        RenderOverrides.SamplingMode() == Sampler.ITERATIVE_SAMPLING) {
      PixelSampler pixels = new TilePixelSampler();
      pixels.setup(x, y, width, height);

      randomSampler = new RandomSampler(x, y, width, height,
                                        sopen, sclose, pixels, 1);
    }
  }

  int maximumSampleCount() {
    return 1;
  }

  int getMoreSamples(List<Sample> samples, RNG rng) {
    if (pass == 0 && randomSampler != null) {
      int count = randomSampler.getMoreSamples(samples, rng);
      if (count != 0) {
        return count;
      }
      pass++;
    }

    while (true) {
      if (currentSample >= wantedSamples) {
        return 0;
      }

      // Generate sample with Halton sequence and reject if outside image extent
      double u = RadicalInverse(currentSample, 3);
      double v = RadicalInverse(currentSample, 2);
      double lerpDelta = Math.max(width, height).toDouble();

      samples[0].imageX = Lerp(u, left, left + lerpDelta);
      samples[0].imageY = Lerp(v, top, top + lerpDelta);

      ++currentSample;

      if (samples[0].imageX > right || samples[0].imageY > bottom) {
        continue;
      }

      break;
    }

    // Generate lens, time, and integrator samples for _HaltonSampler_
    samples[0].lensU = RadicalInverse(currentSample, 5);
    samples[0].lensV = RadicalInverse(currentSample, 7);
    samples[0].time = Lerp(RadicalInverse(currentSample, 11),
                           shutterOpen, shutterClose);

    for (int i = 0; i < samples[0].n1D.length; ++i) {
      LatinHypercube(samples[0].oneD[i], samples[0].n1D[i], 1, rng);
    }

    for (int i = 0; i < samples[0].n2D.length; ++i) {
      LatinHypercube(samples[0].twoD[i], samples[0].n2D[i], 2, rng);
    }

    return 1;
  }

  int roundSize(int size) {
    return size;
  }

  static HaltonSampler Create(ParamSet params, int x, int y, int width,
                              int height, Camera camera, PixelSampler pixels) {
    // Initialize common sampler parameters
    int nsamp = params.findOneInt('pixelsamples', 4);

    return new HaltonSampler(x, y, width, height, camera.shutterOpen,
                             camera.shutterClose, nsamp);
  }

  int wantedSamples;
  int currentSample;
  Sampler randomSampler;
  int pass;
}
