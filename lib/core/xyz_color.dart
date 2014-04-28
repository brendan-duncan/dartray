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

class XYZColor extends Spectrum {
  XYZColor([double v = 0.0])
      : super._(3, v);

  XYZColor.xyz(double x, double y, double z)
      : super._(3) {
    c[0] = x;
    c[1] = y;
    c[2] = z;
  }

  XYZColor.rgb(double r, double g, double b)
      : super._(3) {
    Spectrum.RGBToXYZ(r, g, b, c);
  }

  XYZColor.from(Spectrum color)
      : super._(3) {
    if (color is RGBColor) {
      Spectrum.RGBToXYZ(color.c[0], color.c[1], color.c[2], c);
    } else if (color is XYZColor) {
      c[0] = color.c[0];
      c[1] = color.c[1];
      c[2] = color.c[2];
    } else if (color is SampledSpectrum) {
      Spectrum xyz = color.toXYZ();
      c[0] = xyz.c[0];
      c[1] = xyz.c[1];
      c[2] = xyz.c[2];
    }
  }

  double get x => c[0];

  set x(double v) => c[0] = v;

  double get y => c[1];

  set y(double v) => c[1] = v;

  double get z => c[2];

  set z(double v) => c[2] = v;

  double luminance() => c[1];

  RGBColor toRGB() => new RGBColor.from(this);

  XYZColor toXYZ() => this;

  XYZColor setSampled(List<double> lambda, List<double> v, [int offset = 0]) {
    RGBColor rgb = new RGBColor().setSampled(lambda, v, offset);
    Spectrum.RGBToXYZ(rgb.c[0], rgb.c[1], rgb.c[2], c);
    return this;
  }

  XYZColor setRGB(double r, double g, double b,
                  [int type = Spectrum.SPECTRUM_REFLECTANCE]) {
    Spectrum.RGBToXYZ(r, g, b, c);
    return this;
  }

  XYZColor operator +(XYZColor s) =>
    new XYZColor.xyz(c[0] + s.c[0], c[1] + s.c[1], c[2] + s.c[2]);

  XYZColor operator -(XYZColor s) =>
      new XYZColor.xyz(c[0] - s.c[0], c[1] - s.c[1], c[2] - s.c[2]);

  XYZColor operator *(s) {
    if (s is num) {
      return new XYZColor.xyz(c[0] * s, c[1] * s, c[2] * s);
    }
    if (s is XYZColor) {
      return new XYZColor.xyz(c[0] * s.c[0], c[1] * s.c[1], c[2] * s.c[2]);
    }
    LogSevere('XYZColor or double expected.');
    return new XYZColor(0.0);
  }

  XYZColor operator /(s) {
    if (s is num) {
      return new XYZColor.xyz(c[0] / s, c[1] / s, c[2] / s);
    }
    if (s is XYZColor) {
      return new XYZColor.xyz(c[0] / s.c[0], c[1] / s.c[1], c[2] / s.c[2]);
    }
    LogSevere('XYZColor or double expected.');
    return new XYZColor(0.0);
  }
}
