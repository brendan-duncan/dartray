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

class HaltonSampler extends Sampler {
  HaltonSampler(int xs, int xe, int ys, int ye, double sopen, double sclose,
                int ps, int samplingMode) :
    super(xs, xe, ys, ye, sopen, sclose, ps, samplingMode) {
    int delta = Math.max(xPixelEnd - xPixelStart,
                         yPixelEnd - yPixelStart);
    wantedSamples = samplesPerPixel * delta * delta;
    currentSample = 0;
  }

  int maximumSampleCount() {
    return 1;
  }

  int getMoreSamples(List<Sample> samples, RNG rng) {
    while (true) {
      if (currentSample >= wantedSamples) {
        return 0;
      }

      // Generate sample with Halton sequence and reject if outside image extent
      double u = RadicalInverse(currentSample, 3);
      double v = RadicalInverse(currentSample, 2);
      double lerpDelta = Math.max(xPixelEnd - xPixelStart,
                                  yPixelEnd - yPixelStart).toDouble();

      samples[0].imageX = Lerp(u, xPixelStart, xPixelStart + lerpDelta);
      samples[0].imageY = Lerp(v, yPixelStart, yPixelStart + lerpDelta);

      ++currentSample;

      if (samples[0].imageX >= xPixelEnd || samples[0].imageY >= yPixelEnd) {
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

  Sampler getSubSampler(int num, int count) {
    List<int> range = [0, 0, 0, 0];
    computeSubWindow(num, count, range);
    if (range[0] == range[1] || range[2] == range[3]) {
      return null;
    }
    return new HaltonSampler(range[0], range[1], range[2], range[3],
                             shutterOpen, shutterClose, samplesPerPixel,
                             samplingMode);
  }

  int roundSize(int size) {
    return size;
  }

  static HaltonSampler Create(ParamSet params, Film film, Camera camera,
                              PixelSampler pixels) {
    // Initialize common sampler parameters
    List<int> range = [0, 0, 0, 0];
    film.getSampleExtent(range);
    int nsamp = params.findOneInt('pixelsamples', 4);

    String mode = params.findOneString('mode', 'twopass');
    int samplingMode = (mode == 'full') ? Sampler.FULL_SAMPLING :
                       (mode == 'twopass') ? Sampler.TWO_PASS_SAMPLING :
                       (mode == 'iterative') ? Sampler.ITERATIVE_SAMPLING :
                       -1;
    if (samplingMode == -1) {
      LogWarning('Invalid sampling mode: $mode. Using \'twopass\'.');
      samplingMode = Sampler.TWO_PASS_SAMPLING;
    }

    return new HaltonSampler(range[0], range[1], range[2], range[3],
                             camera.shutterOpen, camera.shutterClose,
                             nsamp, samplingMode);
  }

  int wantedSamples;
  int currentSample;
}
