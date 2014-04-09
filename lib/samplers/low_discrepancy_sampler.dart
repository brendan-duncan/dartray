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
                        int nsamp, double sopen, double sclose) :
    super(xstart, xend, ystart, yend, RoundUpPow2(nsamp), sopen, sclose) {
    //pixels = new LinearImageSampler(xstart, xend, ystart, yend);
    //pixels = new RandomImageSampler(xstart, xend, ystart, yend);
    pixels = new TileImageSampler(xstart, xend, ystart, yend);
    pixelIndex = 0;
    if (!IsPowerOf2(nsamp)) {
      LogWarning('Pixel samples being rounded up to power of 2');
      nPixelSamples = RoundUpPow2(nsamp);
    } else {
      nPixelSamples = nsamp;
    }
    sampleBuf = null;
  }

  static LowDiscrepancySampler Create(ParamSet params, Film film,
                                      Camera camera) {
    // Initialize common sampler parameters
    List<int> extents = [0, 0, 0, 0];
    film.getSampleExtent(extents);
    int nsamp = params.findOneInt('pixelsamples', 4);

    return new LowDiscrepancySampler(extents[0], extents[1], extents[2],
                                     extents[3], nsamp,
                                     camera.shutterOpen, camera.shutterClose);
  }

  Sampler getSubSampler(int num, int count) {
    List extents = [0, 0, 0, 0];
    computeSubWindow(num, count, extents);
    if (extents[0] == extents[1] || extents[2] == extents[3]) {
      return null;
    }

    return new LowDiscrepancySampler(extents[0], extents[1],
                                     extents[2], extents[3],
                                     nPixelSamples,
                                     shutterOpen, shutterClose);
  }

  int roundSize(int size) {
    return RoundUpPow2(size);
  }

  int getMoreSamples(List<Sample> samples, RNG rng) {
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


  ImageSampler pixels;
  Int32List pixel = new Int32List(2);
  int pixelIndex;
  int nPixelSamples;
  Float32List sampleBuf;
}
