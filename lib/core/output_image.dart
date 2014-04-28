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
 * Storage for the output of a renderer.
 *
 * A render thread may work on a portion of the overal image. Only the output
 * of the threads portion of the overal image is stored. The [xOffset],
 * [yOffset], [width], and [height] parameters specify the region of the
 * overal image it occupies.
 *
 * Pixel colors are stored in high dynamic range floating-point values.
 * Tone-mapping will need to be done to translate the pixels into a low dynamic
 * range format.
 */
class OutputImage {
  int width;
  int height;
  int xOffset;
  int yOffset;
  int imageWidth;
  int imageHeight;
  final Float32List rgb;

  OutputImage(this.xOffset, this.yOffset, int width, height,
              [this.imageWidth, this.imageHeight, Float32List rgb])
      : this.width = width,
        this.height = height,
        this.rgb = rgb != null ? rgb : new Float32List(width * height * 3) {
    if (imageWidth == null) {
      imageWidth = width;
    }
    if (imageHeight == null) {
      imageHeight = height;
    }
  }

  Img.Image toImage({double gamma: 2.2}) {
    Img.Image img = new Img.Image(width, height);
    Uint8List pixels = img.getBytes();
    for (int i = 0, oi = 0, len = rgb.length; i < len; i += 3, oi += 4) {
      pixels[oi] = (rgb[i] * 255.0).floor().clamp(0, 255);
      pixels[oi + 1] = (rgb[i + 1] * 255.0).floor().clamp(0, 255);
      pixels[oi + 2] = (rgb[i + 2] * 255.0).floor().clamp(0, 255);
      pixels[oi + 3] = 255;
    }
    return Img.adjustColor(img, gamma: 1.0 / gamma);
  }
}
