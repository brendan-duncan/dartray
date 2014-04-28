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

class SpectrumImage {
  static const int FLOAT = 1;
  static const int SPECTRUM = 3;

  int width;
  int height;
  int samplesPerPixel;
  Float32List data;

  SpectrumImage(int width, int height, [int samplesPerPixel = SPECTRUM])
      : this.width = width,
        this.height = height,
        this.samplesPerPixel = samplesPerPixel,
        data = new Float32List(width * height * samplesPerPixel);

  SpectrumImage.fromImage(Img.Image img)
      : width = img.width,
        height = img.height,
        samplesPerPixel = 3,
        data = new Float32List(img.width * img.height * 3) {
    Uint8List b = img.getBytes();
    for (int y = 0, p = 0, d = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x, p += 4) {
        data[d++] = b[p] / 255.0;
        data[d++] = b[p + 1] / 255.0;
        data[d++] = b[p + 2] / 255.0;
      }
    }
  }

  SpectrumImage.from(SpectrumImage other)
      : width = other.width,
        height = other.height,
        samplesPerPixel = other.samplesPerPixel,
        data = new Float32List.fromList(other.data);

  /**
   * Convert the image to [FLOAT] or [SPECTRUM] format.
   */
  SpectrumImage convert(int format) {
    if (format == samplesPerPixel) {
      return this;
    }

    if (format == FLOAT) {
      SpectrumImage out = new SpectrumImage(width, height, FLOAT);
      int i = 0;
      int j = 0;
      int len = data.length;
      while (i < len) {
        double r = data[i++];
        double g = data[i++];
        double b = data[i++];
        double y = 0.212671 * r + 0.715160 * g + 0.072169 * b;
        out.data[j++] = y;
      }

      return out;
    }

    SpectrumImage out = new SpectrumImage(width, height, SPECTRUM);
    int i = 0;
    int j = 0;
    int len = data.length;
    while (i < len) {
      double y = data[i++];
      out.data[j++] = y;
      out.data[j++] = y;
      out.data[j++] = y;
    }

    return out;
  }

  void set(SpectrumImage other) {
    data = other.data;
    width = other.width;
    height = other.height;
  }

  operator [](int index) {
    if (samplesPerPixel == 1) {
      return data[index];
    }

    index *= samplesPerPixel;
    _output.c[0] = data[index];
    _output.c[1] = data[index + 1];
    _output.c[2] = data[index + 2];
    return _output;
  }

  operator []=(int index, s) {
    if (samplesPerPixel == 1) {
      data[index] = s;
      return;
    }

    index *= samplesPerPixel;
    data[index] = s.c[0];
    data[index + 1] = s.c[1];
    data[index + 2] = s.c[2];
  }

  // Use this to avoid allocating Spectrum for index operator.
  // Dart is non-threaded (other than isolates, which don't share memory).
  static RGBColor _output = new RGBColor(0.0);
}
