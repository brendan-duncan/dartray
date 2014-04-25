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
part of core;

/**
 * Determines the pixels on the film plane for sampling by a [Sampler].
 *
 * This allows the separation of calculating what pixels to render vs what
 * sample within pixels to sample (as done by [Sampler]).
 */
abstract class PixelSampler {
  void setup(int x, int y, int width, int height) {
    this.left = x;
    this.top = y;
    this.width = width;
    this.height = height;
  }

  int get right => left + width - 1;

  int get bottom => top + height - 1;

  int numPixels();

  void getPixel(int index, List<int> pixel);

  int left;
  int top;
  int width;
  int height;
}
