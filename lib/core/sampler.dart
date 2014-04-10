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
part of core;

/**
 * Determines the points on the film plane for tracing rays.
 */
abstract class Sampler {
  Sampler(this.xPixelStart, this.xPixelEnd, this.yPixelStart,
          this.yPixelEnd, this.samplesPerPixel, this.shutterOpen,
          this.shutterClose);

  int get width => (xPixelEnd - xPixelStart);

  int get height => (yPixelEnd - yPixelStart);

  int get top => yPixelStart;

  int get left => xPixelStart;

  int get right => xPixelEnd - 1;

  int get bottom => yPixelEnd - 1;

  int getMoreSamples(List<Sample> sample, RNG rng);

  int maximumSampleCount();

  bool reportResults(List<Sample> samples, List<RayDifferential> rays,
                     List<Spectrum> Ls, List<Intersection> isects,
                     int count) {
    return true;
  }

  Sampler getSubSampler(int num, int count);

  int roundSize(int size);

  void computeSubWindow(int num, int count, List<int> extents) {
    ComputeSubWindow(width, height, num, count, extents);
  }

  static void ComputeSubWindow(int w, int h, int num, int count,
                               List<int> extents) {
    // Determine how many tiles to use in each dimension, nx and ny
    int nx = count;
    int ny = 1;
    while ((nx & 0x1) == 0 && 2 * w * ny < h * nx) {
      nx >>= 1;
      ny <<= 1;
    }
    assert(nx * ny == count);

    // Compute x and y pixel sample range for sub-window
    int xo = num % nx;
    int yo = num ~/ nx;
    double tx0 = xo / nx;
    double tx1 = (xo + 1) / nx;
    double ty0 = yo / ny;
    double ty1 = (yo + 1) / ny;
    extents[0] = Lerp(tx0, 0, w).floor();
    extents[1] = Math.min(Lerp(tx1, 0, w).floor(), w - 1);
    extents[2] = Lerp(ty0, 0, h).floor();
    extents[3] = Math.min(Lerp(ty1, 0, h).floor(), h - 1);
  }

  int xPixelStart;
  int xPixelEnd;
  int yPixelStart;
  int yPixelEnd;
  int samplesPerPixel;
  double shutterOpen;
  double shutterClose;
}

