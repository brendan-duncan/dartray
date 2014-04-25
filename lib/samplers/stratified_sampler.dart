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

class StratifiedSampler extends Sampler {
  StratifiedSampler(int x, int y, int width, int height,
                    this.jitterSamples, double sopen, double sclose,
                    this.pixels, int xs, int ys) :
    super(x, y, width, height, sopen, sclose, xs * ys) {
    if (pixels == null) {
      LogSevere('A PixelSampler is required by StratifiedSampler');
    }
    pixels.setup(x, y, width, height);
    pixelIndex = 0;
    xPixelSamples = xs;
    yPixelSamples = ys;
    nPixelSamples = xPixelSamples * yPixelSamples;
    imageSamples = new Float32List(2 * nPixelSamples);
    lensSamples = new Float32List(2 * nPixelSamples);
    timeSamples = new Float32List(xPixelSamples * yPixelSamples);

    pass = 0;
    if (RenderOverrides.SamplingMode() == Sampler.TWO_PASS_SAMPLING ||
        RenderOverrides.SamplingMode() == Sampler.ITERATIVE_SAMPLING) {
      randomSampler = new RandomSampler(x, y, width, height,
                                        sopen, sclose, pixels, 1);
    }
  }

  int roundSize(int size) {
    return size;
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

    // Generate stratified camera samples

    // Generate initial stratified samples into _sampleBuf_ memory
    int nb = 2 * nPixelSamples;

    StratifiedSample2D(imageSamples, xPixelSamples, yPixelSamples, rng,
                       jitterSamples);

    StratifiedSample2D(lensSamples, xPixelSamples, yPixelSamples, rng,
                       jitterSamples);

    StratifiedSample1D(timeSamples, xPixelSamples * yPixelSamples, rng,
                       jitterSamples);

    pixels.getPixel(pixelIndex++, pixel);

    // Shift stratified image samples to pixel coordinates
    for (int o = 0; o < 2 * xPixelSamples * yPixelSamples; o += 2) {
      imageSamples[o] += pixel[0];
      imageSamples[o + 1] += pixel[1];
    }

    // Decorrelate sample dimensions
    Shuffle(lensSamples, 0, xPixelSamples * yPixelSamples, 2, rng);
    Shuffle(timeSamples, 0, xPixelSamples * yPixelSamples, 1, rng);

    // Initialize stratified _samples_ with sample values
    for (int i = 0; i < nPixelSamples; ++i) {
      samples[i].imageX = imageSamples[2 * i];
      samples[i].imageY = imageSamples[2 * i + 1];
      samples[i].lensU = lensSamples[2 * i];
      samples[i].lensV = lensSamples[2 * i + 1];
      samples[i].time = Lerp(timeSamples[i], shutterOpen, shutterClose);
      // Generate stratified samples for integrators
      for (int j = 0; j < samples[i].n1D.length; ++j) {
        LatinHypercube(samples[i].oneD[j], samples[i].n1D[j], 1, rng);
      }

      for (int j = 0; j < samples[i].n2D.length; ++j) {
        LatinHypercube(samples[i].twoD[j], samples[i].n2D[j], 2, rng);
      }
    }

    return nPixelSamples;
  }

  int maximumSampleCount() {
    return nPixelSamples;
  }

  static StratifiedSampler Create(ParamSet params, int x, int y, int width,
                                  int height, Camera camera,
                                  PixelSampler pixels) {
    bool jitter = params.findOneBool('jitter', true);

    int pixelsamples = params.findOneInt('pixelsamples', null);

    int xsamp;
    int ysamp;
    if (pixelsamples != null) {
      xsamp = pixelsamples;
      ysamp = pixelsamples;
    } else {
      xsamp = params.findOneInt('xsamples', 2);
      ysamp = params.findOneInt('ysamples', 2);
    }

    return new StratifiedSampler(x, y, width, height, jitter,
                                 camera.shutterOpen, camera.shutterClose,
                                 pixels, xsamp, ysamp);
  }

  int xPixelSamples;
  int yPixelSamples;
  int nPixelSamples;
  bool jitterSamples;
  PixelSampler pixels;
  Int32List pixel = new Int32List(2);
  int pixelIndex;
  Float32List imageSamples;
  Float32List lensSamples;
  Float32List timeSamples;
  Sampler randomSampler;
  int pass;
}
