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

class SampledSpectrum extends Spectrum {
  static const int NUM_SAMPLES = 4;

  SampledSpectrum([double v = 0.0])
      : super._(NUM_SAMPLES, v);

  SampledSpectrum.from(Spectrum s, [int type = Spectrum.SPECTRUM_REFLECTANCE])
      : super._(NUM_SAMPLES) {
    if (s is SampledSpectrum) {
      for (int i = 0; i < NUM_SAMPLES; ++i) {
        c[i] = s.c[i];
      }
    } else if (s is RGBColor) {
      setRGB(s.c[0], s.c[1], s.c[2], type);
    } else if (s is XYZColor) {
      RGBColor rgb = new RGBColor.from(s);
      setRGB(rgb.c[0], rgb.c[1], rgb.c[2], type);
    }
  }

  SampledSpectrum.rgb(double r, double g, double b)
      : super._(NUM_SAMPLES) {
    setRGB(r, g, b);
  }

  SampledSpectrum.xyz(double x, double y, double z)
      : super._(NUM_SAMPLES) {
    setXYZ(x, y, z);
  }

  SampledSpectrum.fromSampled(List<double> lambda, List<double> v,
                             [int offset = 0])
      : super._(NUM_SAMPLES) {
    setSampled(lambda, v, offset);
  }

  SampledSpectrum setSampled(List<double> lambda, List<double> v,
                             [int offset = 0]) {
    // Sort samples if unordered, use sorted for returned spectrum
    if (!Spectrum.SpectrumSamplesSorted(lambda)) {
      Spectrum.SortSpectrumSamples(lambda, v, offset);
    }

    for (int i = 0; i < NUM_SAMPLES; ++i) {
      // Compute average value of given SPD over $i$th sample's range
      double lambda0 = Lerp(i / NUM_SAMPLES,
                            Spectrum.SAMPLED_LAMBDA_START,
                            Spectrum.SAMPLED_LAMBDA_END);

      double lambda1 = Lerp((i + 1) / NUM_SAMPLES,
          Spectrum.SAMPLED_LAMBDA_START, Spectrum.SAMPLED_LAMBDA_END);

      c[i] = Spectrum.AverageSpectrumSamples(lambda, v, lambda0, lambda1,
                                             offset);
    }
    return this;
  }

  SampledSpectrum operator +(SampledSpectrum s) {
    SampledSpectrum r = new SampledSpectrum();
    for (int i = 0; i < NUM_SAMPLES; ++i) {
      r.c[i] = c[i] + s.c[i];
    }
    return r;
  }

  SampledSpectrum operator -(SampledSpectrum s) {
    SampledSpectrum r = new SampledSpectrum();
    for (int i = 0; i < NUM_SAMPLES; ++i) {
      r.c[i] = c[i] - s.c[i];
    }
    return r;
  }

  SampledSpectrum operator *(s) {
    if (s is num) {
      SampledSpectrum r = new SampledSpectrum();
      for (int i = 0; i < NUM_SAMPLES; ++i) {
        r.c[i] = c[i] * s;
      }
      return r;
    }
    if (s is SampledSpectrum) {
      SampledSpectrum r = new SampledSpectrum();
      for (int i = 0; i < NUM_SAMPLES; ++i) {
        r.c[i] = c[i] * s.c[i];
      }
      return r;
    }
    LogSevere('SampledSpectrum or num expected.');
    return new SampledSpectrum(0.0);
  }

  SampledSpectrum operator /(s) {
    if (s is num) {
      SampledSpectrum r = new SampledSpectrum();
      for (int i = 0; i < NUM_SAMPLES; ++i) {
        r.c[i] = c[i] / s;
      }
      return r;
    }
    if (s is SampledSpectrum) {
      SampledSpectrum r = new SampledSpectrum();
      for (int i = 0; i < NUM_SAMPLES; ++i) {
        r.c[i] = c[i] / s.c[i];
      }
      return r;
    }
    LogSevere('SampledSpectrum or double expected.');
    return new SampledSpectrum(0.0);
  }

  SampledSpectrum setXYZ(double x, double y, double z) {
    RGBColor rgb = new RGBColor.xyz(x, y, z);
    setRGB(rgb.c[0], rgb.c[1], rgb.c[2]);
    return this;
  }

  SampledSpectrum setRGB(double r, double g, double b,
              [int type = Spectrum.SPECTRUM_REFLECTANCE]) {
    SampledSpectrum res = new SampledSpectrum();

    if (type == Spectrum.SPECTRUM_REFLECTANCE) {
      // Convert reflectance spectrum to RGB
      if (r <= g && r <= b) {
        // Compute reflectance _SampledSpectrum_ with _r_ as minimum
        res += _Spectrum.G.rgbRefl2SpectWhite * r;
        if (g <= b) {
          res += _Spectrum.G.rgbRefl2SpectCyan * (g - r);
          res += _Spectrum.G.rgbRefl2SpectBlue * (b - g);
        } else {
          res += _Spectrum.G. rgbRefl2SpectCyan * (b - r);
          res += _Spectrum.G.rgbRefl2SpectGreen * (g - b);
        }
      } else if (g <= r && g <= b) {
        // Compute reflectance _SampledSpectrum_ with _g_ as minimum
        res += _Spectrum.G.rgbRefl2SpectWhite * g;
        if (r <= b) {
          res += _Spectrum.G.rgbRefl2SpectMagenta * (r - g);
          res += _Spectrum.G.rgbRefl2SpectBlue * (b - r);
        } else {
          res += _Spectrum.G.rgbRefl2SpectMagenta * (b - g);
          res += _Spectrum.G.rgbRefl2SpectRed * (r - b);
        }
      } else {
        // Compute reflectance _SampledSpectrum_ with _b_ as minimum
        res += _Spectrum.G.rgbRefl2SpectWhite * b;
        if (r <= g) {
          res += _Spectrum.G.rgbRefl2SpectYellow * (r - b);
          res += _Spectrum.G.rgbRefl2SpectGreen * (g - r);
        } else {
          res += _Spectrum.G.rgbRefl2SpectYellow * (g - b);
          res += _Spectrum.G.rgbRefl2SpectRed * (r - g);
        }
      }

      res *= 0.94;
    } else {
      // Convert illuminant spectrum to RGB
      if (r <= g && r <= b) {
        // Compute illuminant _SampledSpectrum_ with _r_ as minimum
        res += _Spectrum.G.rgbIllum2SpectWhite * r;
        if (g <= b) {
          res += _Spectrum.G.rgbIllum2SpectCyan * (g - r);
          res += _Spectrum.G.rgbIllum2SpectBlue * (b - g);
        } else {
          res += _Spectrum.G.rgbIllum2SpectCyan * (b - r);
          res += _Spectrum.G.rgbIllum2SpectGreen * (g - b);
        }
      } else if (g <= r && g <= b) {
        // Compute illuminant _SampledSpectrum_ with _g_ as minimum
        res += _Spectrum.G.rgbIllum2SpectWhite * g;
        if (r <= b) {
          res += _Spectrum.G.rgbIllum2SpectMagenta * (r - g);
          res += _Spectrum.G.rgbIllum2SpectBlue * (b - r);
        } else {
          res += _Spectrum.G.rgbIllum2SpectMagenta * (b - g);
          res += _Spectrum.G.rgbIllum2SpectRed * (r - b);
        }
      } else {
        // Compute illuminant _SampledSpectrum_ with _b_ as minimum
        res += _Spectrum.G.rgbIllum2SpectWhite * b;
        if (r <= g) {
          res += _Spectrum.G.rgbIllum2SpectYellow * (r - b);
          res += _Spectrum.G.rgbIllum2SpectGreen * (g - r);
        } else {
          res += _Spectrum.G.rgbIllum2SpectYellow * (g - b);
          res += _Spectrum.G.rgbIllum2SpectRed * (r - g);
        }
      }

      res *= 0.86445;
    }

    c[0] = res.c[0].clamp(0.0, INFINITY);
    c[1] = res.c[1].clamp(0.0, INFINITY);
    c[2] = res.c[2].clamp(0.0, INFINITY);

    return this;
  }

  XYZColor toXYZ() {
    XYZColor xyz = new XYZColor();
    xyz.c[0] = 0.0;
    xyz.c[1] = 0.0;
    xyz.c[2] = 0.0;

    for (int i = 0; i < NUM_SAMPLES; ++i) {
      xyz.c[0] += _Spectrum.G.X.c[i] * c[i];
      xyz.c[1] += _Spectrum.G.Y.c[i] * c[i];
      xyz.c[2] += _Spectrum.G.Z.c[i] * c[i];
    }

    double scale = (Spectrum.SAMPLED_LAMBDA_END -
                    Spectrum.SAMPLED_LAMBDA_START) /
                   (Spectrum.CIE_Y_integral * NUM_SAMPLES);

    xyz.c[0] *= scale;
    xyz.c[1] *= scale;
    xyz.c[2] *= scale;

    return xyz;
  }

  double luminance() {
    double yy = 0.0;
    for (int i = 0; i < NUM_SAMPLES; ++i) {
      yy += _Spectrum.G.Y.c[i] * c[i];
    }

    return yy * (Spectrum.SAMPLED_LAMBDA_END -
                 Spectrum.SAMPLED_LAMBDA_START) /
                (Spectrum.CIE_Y_integral * NUM_SAMPLES);
  }

  RGBColor toRGB() {
    XYZColor xyz = toXYZ();
    return new RGBColor.from(xyz);
  }

  SampledSpectrum clamp([double low = 0.0, double high = INFINITY]) {
    SampledSpectrum r = new SampledSpectrum();
    for (int i = 0; i < NUM_SAMPLES; ++i) {
      r.c[i] = c[i].clamp(low, high);
    }
    return r;
  }
}

class _Spectrum {
  static _Spectrum G = new _Spectrum();

  SampledSpectrum X, Y, Z;
  SampledSpectrum rgbRefl2SpectWhite, rgbRefl2SpectCyan;
  SampledSpectrum rgbRefl2SpectMagenta, rgbRefl2SpectYellow;
  SampledSpectrum rgbRefl2SpectRed, rgbRefl2SpectGreen;
  SampledSpectrum rgbRefl2SpectBlue;
  SampledSpectrum rgbIllum2SpectWhite, rgbIllum2SpectCyan;
  SampledSpectrum rgbIllum2SpectMagenta, rgbIllum2SpectYellow;
  SampledSpectrum rgbIllum2SpectRed, rgbIllum2SpectGreen;
  SampledSpectrum rgbIllum2SpectBlue;

  _Spectrum() {
    X = new SampledSpectrum();
    Y = new SampledSpectrum();
    Z = new SampledSpectrum();
    rgbRefl2SpectWhite = new SampledSpectrum();
    rgbRefl2SpectCyan = new SampledSpectrum();
    rgbRefl2SpectMagenta = new SampledSpectrum();
    rgbRefl2SpectYellow = new SampledSpectrum();
    rgbRefl2SpectRed = new SampledSpectrum();
    rgbRefl2SpectGreen = new SampledSpectrum();
    rgbRefl2SpectBlue = new SampledSpectrum();
    rgbIllum2SpectWhite = new SampledSpectrum();
    rgbIllum2SpectCyan = new SampledSpectrum();
    rgbIllum2SpectMagenta = new SampledSpectrum();
    rgbIllum2SpectYellow = new SampledSpectrum();
    rgbIllum2SpectRed = new SampledSpectrum();
    rgbIllum2SpectGreen = new SampledSpectrum();
    rgbIllum2SpectBlue = new SampledSpectrum();

    // Compute XYZ matching functions for [SampledSpectrum]
    for (int i = 0; i < SampledSpectrum.NUM_SAMPLES; ++i) {
      double wl0 = Lerp(i / SampledSpectrum.NUM_SAMPLES,
          Spectrum.SAMPLED_LAMBDA_START, Spectrum.SAMPLED_LAMBDA_END);

      double wl1 = Lerp((i + 1) / SampledSpectrum.NUM_SAMPLES,
          Spectrum.SAMPLED_LAMBDA_START, Spectrum.SAMPLED_LAMBDA_END);

      X.c[i] = Spectrum.AverageSpectrumSamples(Spectrum.CIE_lambda,
                                               Spectrum.CIE_X,
                                               wl0, wl1);
      Y.c[i] = Spectrum.AverageSpectrumSamples(Spectrum.CIE_lambda,
                                               Spectrum.CIE_Y,
                                               wl0, wl1);
      Z.c[i] = Spectrum.AverageSpectrumSamples(Spectrum.CIE_lambda,
                                               Spectrum.CIE_Z,
                                               wl0, wl1);
    }

    // Compute RGB to spectrum functions for _SampledSpectrum_
    for (int i = 0; i < SampledSpectrum.NUM_SAMPLES; ++i) {
        double wl0 = Lerp(i / SampledSpectrum.NUM_SAMPLES,
            Spectrum.SAMPLED_LAMBDA_START, Spectrum.SAMPLED_LAMBDA_END);

        double wl1 = Lerp((i + 1) / SampledSpectrum.NUM_SAMPLES,
            Spectrum.SAMPLED_LAMBDA_START, Spectrum.SAMPLED_LAMBDA_END);

        rgbRefl2SpectWhite.c[i] =
            Spectrum.AverageSpectrumSamples(Spectrum.RGB2SpectLambda,
                                            Spectrum.RGBRefl2SpectWhite,
                                            wl0, wl1);
        rgbRefl2SpectCyan.c[i] =
            Spectrum.AverageSpectrumSamples(Spectrum.RGB2SpectLambda,
                                            Spectrum.RGBRefl2SpectCyan,
                                            wl0, wl1);
        rgbRefl2SpectMagenta.c[i] =
            Spectrum.AverageSpectrumSamples(Spectrum.RGB2SpectLambda,
                                            Spectrum.RGBRefl2SpectMagenta,
                                            wl0, wl1);
        rgbRefl2SpectYellow.c[i] =
            Spectrum.AverageSpectrumSamples(Spectrum.RGB2SpectLambda,
                                            Spectrum.RGBRefl2SpectYellow,
                                            wl0, wl1);
        rgbRefl2SpectRed.c[i] =
            Spectrum.AverageSpectrumSamples(Spectrum.RGB2SpectLambda,
                                            Spectrum.RGBRefl2SpectRed,
                                            wl0, wl1);
        rgbRefl2SpectGreen.c[i] =
            Spectrum.AverageSpectrumSamples(Spectrum.RGB2SpectLambda,
                                            Spectrum.RGBRefl2SpectGreen,
                                            wl0, wl1);
        rgbRefl2SpectBlue.c[i] =
            Spectrum.AverageSpectrumSamples(Spectrum.RGB2SpectLambda,
                                            Spectrum.RGBRefl2SpectBlue,
                                            wl0, wl1);

        rgbIllum2SpectWhite.c[i] =
            Spectrum.AverageSpectrumSamples(Spectrum.RGB2SpectLambda,
                                            Spectrum.RGBIllum2SpectWhite,
                                            wl0, wl1);
        rgbIllum2SpectCyan.c[i] =
            Spectrum.AverageSpectrumSamples(Spectrum.RGB2SpectLambda,
                                            Spectrum.RGBIllum2SpectCyan,
                                            wl0, wl1);
        rgbIllum2SpectMagenta.c[i] =
            Spectrum.AverageSpectrumSamples(Spectrum.RGB2SpectLambda,
                                            Spectrum.RGBIllum2SpectMagenta,
                                            wl0, wl1);
        rgbIllum2SpectYellow.c[i] =
            Spectrum.AverageSpectrumSamples(Spectrum.RGB2SpectLambda,
                                            Spectrum.RGBIllum2SpectYellow,
                                            wl0, wl1);
        rgbIllum2SpectRed.c[i] =
            Spectrum.AverageSpectrumSamples(Spectrum.RGB2SpectLambda,
                                            Spectrum.RGBIllum2SpectRed,
                                            wl0, wl1);
        rgbIllum2SpectGreen.c[i] =
            Spectrum.AverageSpectrumSamples(Spectrum.RGB2SpectLambda,
                                            Spectrum.RGBIllum2SpectGreen,
                                            wl0, wl1);
        rgbIllum2SpectBlue.c[i] =
            Spectrum.AverageSpectrumSamples(Spectrum.RGB2SpectLambda,
                                            Spectrum.RGBIllum2SpectBlue,
                                            wl0, wl1);
    }
  }
}
