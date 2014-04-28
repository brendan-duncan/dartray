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

class RGBColor extends Spectrum {
  RGBColor([double v = 0.0])
      : super._(3, v);

  RGBColor.rgb(double r, double g, double b)
      : super._(3) {
    c[0] = r;
    c[1] = g;
    c[2] = b;
    assert(!c[0].isNaN && !c[1].isNaN && !c[2].isNaN);
    assert(!c[0].isInfinite && !c[1].isInfinite && !c[2].isInfinite);
  }

  RGBColor.xyz(double x, double y, double z)
      : super._(3) {
    Spectrum.XYZToRGB(x, y, z, c);
    assert(!c[0].isNaN && !c[1].isNaN && !c[2].isNaN);
    assert(!c[0].isInfinite && !c[1].isInfinite && !c[2].isInfinite);
  }

  RGBColor.from(Spectrum s)
      : super._(3) {
    if (s is RGBColor) {
      c[0] = s.c[0];
      c[1] = s.c[1];
      c[2] = s.c[2];
    } else if (s is XYZColor) {
      Spectrum.XYZToRGB(s.c[0], s.c[1], s.c[2], c);
    } else if (s is SampledSpectrum) {
      Spectrum rgb = s.toRGB();
      c[0] = rgb.c[0];
      c[1] = rgb.c[1];
      c[2] = rgb.c[2];
    }
    assert(!c[0].isNaN && !c[1].isNaN && !c[2].isNaN);
    assert(!c[0].isInfinite && !c[1].isInfinite && !c[2].isInfinite);
  }

  RGBColor setSampled(List<double> lambda, List<double> v, [int offset = 0]) {
    // Sort samples if unordered, use sorted for returned spectrum
    if (!Spectrum.SpectrumSamplesSorted(lambda)) {
      Spectrum.SortSpectrumSamples(lambda, v, offset);
    }

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double yint = 0.0;

    for (int i = 0; i < Spectrum.NUM_CIE_SAMPLES; ++i) {
      yint += Spectrum.CIE_Y[i];
      double val = Spectrum.InterpolateSpectrumSamples(lambda, v,
                                                       Spectrum.CIE_lambda[i],
                                                       offset);
      x += (val * Spectrum.CIE_X[i]);
      y += (val * Spectrum.CIE_Y[i]);
      z += (val * Spectrum.CIE_Z[i]);
    }

    x /= yint;
    y /= yint;
    z /= yint;

    Spectrum.XYZToRGB(x, y, z, c);

    assert(!c[0].isNaN && !c[1].isNaN && !c[2].isNaN);
    assert(!c[0].isInfinite && !c[1].isInfinite && !c[2].isInfinite);

    return this;
  }

  double get r => c[0];

  set r(double v) => c[0] = v;

  double get g => c[1];

  set g(double v) => c[1] = v;

  double get b => c[2];

  set b(double v) => c[2] = v;

  void setXYZ(double x, double y, double z) {
    Spectrum.XYZToRGB(x, y, z, c);
    assert(!c[0].isNaN && !c[1].isNaN && !c[2].isNaN);
    assert(!c[0].isInfinite && !c[1].isInfinite && !c[2].isInfinite);
  }

  XYZColor toXYZ() {
    return new XYZColor.from(this);
  }

  RGBColor toRGB() => this;

  void set(double v) {
    c[0] = v;
    c[1] = v;
    c[2] = v;
    assert(!c[0].isNaN && !c[1].isNaN && !c[2].isNaN);
    assert(!c[0].isInfinite && !c[1].isInfinite && !c[2].isInfinite);
  }

  RGBColor setRGB(double r, double g, double b,
                  [int type = Spectrum.SPECTRUM_REFLECTANCE]) {
    c[0] = r;
    c[1] = g;
    c[2] = b;
    assert(!c[0].isNaN && !c[1].isNaN && !c[2].isNaN);
    assert(!c[0].isInfinite && !c[1].isInfinite && !c[2].isInfinite);
    return this;
  }

  RGBColor operator+(RGBColor s) =>
    new RGBColor.rgb(c[0] + s.c[0], c[1] + s.c[1], c[2] + s.c[2]);

  RGBColor operator-(RGBColor s) =>
      new RGBColor.rgb(c[0] - s.c[0], c[1] - s.c[1], c[2] - s.c[2]);

  RGBColor operator*(s) {
    if (s is num) {
      return new RGBColor.rgb(c[0] * s, c[1] * s, c[2] * s);
    }
    if (s is RGBColor) {
      return new RGBColor.rgb(c[0] * s.c[0], c[1] * s.c[1], c[2] * s.c[2]);
    }
    LogSevere('RGBSpectrum or double expected.');
    return new RGBColor(0.0);
  }

  RGBColor operator/(s) {
    if (s is num) {
      return new RGBColor.rgb(c[0] / s, c[1] / s, c[2] / s);
    }
    if (s is RGBColor) {
      return new RGBColor.rgb(c[0] / s.c[0], c[1] / s.c[1], c[2] / s.c[2]);
    }
    LogSevere('RGBSpectrum or double expected.');
    return new RGBColor(0.0);
  }

  RGBColor operator-() =>
        new RGBColor.rgb(-c[0], -c[1], -c[2]);

  double luminance() {
    return 0.212671 * c[0] + 0.715160 * c[1] + 0.072169 * c[2];
  }

  bool isBlack() {
    if (c[0] != 0.0 || c[1] != 0.0 || c[2] != 0.0) {
      return false;
    }
    return true;
  }

  RGBColor sqrt() =>
    new RGBColor.rgb(Math.sqrt(c[0]), Math.sqrt(c[1]), Math.sqrt(c[2]));

  RGBColor pow(double e) =>
      new RGBColor.rgb(Math.pow(c[0], e),
                       Math.pow(c[1], e),
                       Math.pow(c[2], e));

  RGBColor exp() =>
      new RGBColor.rgb(Math.exp(c[0]), Math.exp(c[1]), Math.exp(c[2]));

  RGBColor clamp([double low = 0.0, double high = INFINITY]) =>
      new RGBColor.rgb(c[0].clamp(low, high),
                       c[1].clamp(low, high),
                       c[2].clamp(low, high));

  String toString() {
    return '${c[0]} ${c[1]} ${c[2]}';
  }

  static RGBColor Lerp(double t, RGBColor s1, RGBColor s2) {
    return (s1 * (1.0 - t)) + (s2 * t);
  }

  static List<RGBColor> AllocateList(int length) {
    List<RGBColor> l = new List<RGBColor>(length);
    for (int i = 0; i < length; ++i) {
      l[i] = new RGBColor(0.0);
    }
    return l;
  }
}
