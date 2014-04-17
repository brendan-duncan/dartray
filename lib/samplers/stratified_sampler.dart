/****************************************************************************
 *  Copyright (C) 2014 by Brendan Duncan.                                   *
 *                                                                          *
 *  This file is part of DartRay.                                           *
 *                                                                          *
 *  Licensed under the Apache License, Version 2.0 (the 'License');         *
 *  you may not use this file except in compliance with the License.        *
 *  You may obtain a copy of the License at                                 *
 *                                                                          *
 *  http://www.apache.org/licenses/LICENSE-2.0                              *
 *                                                                          *
 *  Unless required by applicable law or agreed to in writing, software     *
 *  distributed under the License is distributed on an 'AS IS' BASIS,       *
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 *  See the License for the specific language governing permissions and     *
 *  limitations under the License.                                          *
 *                                                                          *
 *   This project is based on PBRT v2 ; see http://www.pbrt.org             *
 *   pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.*
 ****************************************************************************/
part of samplers;

class StratifiedSampler extends Sampler {
  static StratifiedSampler Create(ParamSet params, Film film, Camera camera,
                                  PixelSampler pixels) {
    bool jitter = params.findOneBool('jitter', true);
    // Initialize common sampler parameters
    List<int> extents = [0, 0, 0, 0];
    film.getSampleExtent(extents);

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

    return new StratifiedSampler(extents[0], extents[1], extents[2],
                                 extents[3], xsamp, ysamp, jitter,
                                 camera.shutterOpen, camera.shutterClose,
                                 pixels);
  }

  StratifiedSampler(int xstart, int xend, int ystart, int yend,
                    int xs, int ys, this.jitterSamples,
                    double sopen, double sclose,
                    this.pixels) :
    super(xstart, xend, ystart, yend, xs * ys, sopen, sclose) {
    if (pixels == null) {
      LogSevere('Pixel sampler is required by StratifiedSampler');
    }
    pixels.setup(xstart, xend, ystart, yend);
    pixelIndex = 0;
    xPixelSamples = xs;
    yPixelSamples = ys;
    nPixelSamples = xPixelSamples * yPixelSamples;
    imageSamples = new Float32List(2 * nPixelSamples);
    lensSamples = new Float32List(2 * nPixelSamples);
    timeSamples = new Float32List(xPixelSamples * yPixelSamples);
  }

  int roundSize(int size) {
    return size;
  }

  Sampler getSubSampler(int num, int count) {
    List<int> range = [0, 0, 0, 0];
    computeSubWindow(num, count, range);
    if (range[0] == range[1] || range[2] == range[3]) {
      return null;
    }
    return new StratifiedSampler(range[0], range[1], range[2], range[3],
                                 xPixelSamples, yPixelSamples, jitterSamples,
                                 shutterOpen, shutterClose, pixels);
  }

  int getMoreSamples(List<Sample> samples, RNG rng) {
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
}
