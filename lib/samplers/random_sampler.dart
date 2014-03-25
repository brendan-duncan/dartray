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
  static RandomSampler Create(ParamSet params, Film film, Camera camera) {
    int ns = params.findOneInt("pixelsamples", 4);
    List<int> extents = [0, 0, 0, 0];
    film.getSampleExtent(extents);
    return new RandomSampler(extents[0], extents[1], extents[2], extents[3], ns,
                             camera.shutterOpen, camera.shutterClose);
  }

  RandomSampler(int xstart, int xend, int ystart,
      int yend, int ns, double sopen, double sclose) :
    super(xstart, xend, ystart, yend, ns, sopen, sclose) {
    xPos = xPixelStart;
    yPos = yPixelStart;
    nSamples = ns;
    // Get storage for a pixel's worth of stratified samples
    imageSamples = new Float32List(2 * nSamples);
    lensSamples = new Float32List(2 * nSamples);
    timeSamples = new Float32List(nSamples);

    RNG rng = new RNG(xstart + ystart * (xend - xstart));
    int i = 0;
    for (; i < nSamples; ++i) {
      timeSamples[i] = rng.randomFloat();
    }
    for (; i < 2 * nSamples; ++i) {
      imageSamples[i] = rng.randomFloat();
      lensSamples[i] = rng.randomFloat();
    }

    // Shift image samples to pixel coordinates
    for (int o = 0, len = 2 * nSamples; o < len; o += 2) {
      imageSamples[o] += xPos;
      imageSamples[o+1] += yPos;
    }

    samplePos = 0;
  }

  int maximumSampleCount() {
    return 1;
  }

  int getMoreSamples(List<Sample> sample, RNG rng) {
    if (samplePos == nSamples) {
      if (xPixelStart == xPixelEnd || yPixelStart == yPixelEnd) {
        return 0;
      }

      if (++xPos == xPixelEnd) {
        xPos = xPixelStart;
        ++yPos;
      }

      if (yPos == yPixelEnd) {
        return 0;
      }

      int i = 0;
      for (; i < nSamples; ++i) {
        timeSamples[i] = rng.randomFloat();
      }
      for (; i < 2 * nSamples; ++i) {
        imageSamples[i] = rng.randomFloat();
        lensSamples[i] = rng.randomFloat();
      }

      // Shift image samples to pixel coordinates
      for (int o = 0; o < 2 * nSamples; o += 2) {
        imageSamples[o] += xPos;
        imageSamples[o + 1] += yPos;
      }

      samplePos = 0;
    }

    // Return next sample point
    sample[0].imageX = imageSamples[2 * samplePos];
    sample[0].imageY = imageSamples[2 * samplePos + 1];
    sample[0].lensU = lensSamples[2 * samplePos];
    sample[0].lensV = lensSamples[2 * samplePos + 1];
    sample[0].time = Lerp(timeSamples[samplePos], shutterOpen, shutterClose);

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

    ++samplePos;
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
                             nSamples, shutterOpen, shutterClose);
  }

  int xPos, yPos, nSamples;
  Float32List imageSamples;
  Float32List lensSamples;
  Float32List timeSamples;
  int samplePos;
}
