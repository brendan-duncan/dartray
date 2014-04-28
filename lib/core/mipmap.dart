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

class MIPMap {
  static const int TEXTURE_REPEAT = 0;
  static const int TEXTURE_BLACK = 1;
  static const int TEXTURE_CLAMP = 2;

  static String GetTextureName(String filename,
                               {bool doTri: false, double maxAniso: 8.0,
                                int wrap: TEXTURE_REPEAT,
                                scale: 1.0,
                                double gamma: 1.0,
                                bool spectrum: true}) {
    String name = filename;
    if (doTri) {
      name += '_TRI:$doTri';
    }
    if (maxAniso != 8.0) {
      name += '_ANI:$maxAniso';
    }
    if (wrap != TEXTURE_REPEAT) {
      name += '_WRAP:$wrap';
    }
    if (scale is num && scale != 1.0) {
      name += '_SCALE:$scale';
    }
    if (scale is Spectrum && !scale.isValue(1.0)) {
      name += '_SCALE:$scale';
    }
    if (gamma != 1.0) {
      name += '_GAMMA:$gamma';
    }
    if (!spectrum) {
      name += '_SPECTRUM:$spectrum';
    }
    return name;
  }

  MIPMap()
      : width = 0,
        height = 0,
        levels = 0;

  MIPMap.texture(SpectrumImage img, String filename,
                 [this.doTrilinear = false,
                  this.maxAnisotropy = 8.0,
                  this.wrapMode = TEXTURE_REPEAT]) {
    int xres = img.width;
    int yres = img.height;

    SpectrumImage resampledImage;
    if (!IsPowerOf2(xres) || !IsPowerOf2(yres)) {
      // Resample image to power-of-two resolution
      int sPow2 = RoundUpPow2(xres);
      int tPow2 = RoundUpPow2(yres);
      if (filename.isNotEmpty) {
        LogInfo('Resizing Image $filename to $sPow2 $tPow2');
      }

      // Resample image in s direction
      List<_ResampleWeight> sWeights = _resampleWeights(xres, sPow2);
      resampledImage = new SpectrumImage(sPow2, tPow2, img.samplesPerPixel);

      var zero = img.samplesPerPixel == 1 ? 0.0 : new Spectrum(0.0);

      // Apply _sWeights_ to zoom in s direction
      for (int t = 0, p = 0; t < yres; ++t) {
        for (int s = 0; s < sPow2; ++s, ++p) {
          // Compute texel (s,t) in s-zoomed image
          resampledImage[p] = zero;
          for (int j = 0; j < 4; ++j) {
            int origS = sWeights[s].firstTexel + j;
            if (wrapMode == TEXTURE_REPEAT) {
              origS = origS % xres;
            } else if (wrapMode == TEXTURE_CLAMP) {
              origS = origS.clamp(0, xres - 1);
            }

            if (origS >= 0 && origS < xres) {
              var px = img[t * xres + origS] *
                       sWeights[s].weight[j];
              resampledImage[t * sPow2 + s] += px;
            }
          }
        }
      }

      // Resample image in $t$ direction
      List<_ResampleWeight> tWeights = _resampleWeights(yres, tPow2);
      List workData = new List(tPow2);

      for (int s = 0; s < sPow2; ++s) {
        for (int t = 0; t < tPow2; ++t) {
          workData[t] = img.samplesPerPixel == 3 ? new Spectrum(0.0) : 0.0;
          for (int j = 0; j < 4; ++j) {
            int offset = tWeights[t].firstTexel + j;
            if (wrapMode == TEXTURE_REPEAT) {
              offset = offset % yres;
            } else if (wrapMode == TEXTURE_CLAMP) {
              offset = offset.clamp(0, yres - 1);
            }

            if (offset >= 0 && offset < yres) {
              var px = resampledImage[offset * sPow2 + s] *
                       tWeights[t].weight[j];
              workData[t] += px;
            }
          }
        }

        for (int t = 0; t < tPow2; ++t) {
          resampledImage[t * sPow2 + s] = _clamp(workData[t]);
        }
      }

      img.set(resampledImage);
      xres = sPow2;
      yres = tPow2;
    }

    width = xres;
    height = yres;
    // Initialize levels of MIPMap from image
    levels = 1 + Log2(Math.max(xres, yres)).toInt();
    pyramid = new List<SpectrumImage>(levels);
    if (filename.isNotEmpty && levels > 1) {
      LogInfo('$filename: Generating $levels MIPMap Levels');
    }

    // Initialize most detailed level of MIPMap
    pyramid[0] = new SpectrumImage.from(img);

    for (int i = 1; i < levels; ++i) {
      // Initialize i'th MIPMap level from i-1'th level
      int sRes = Math.max(1, pyramid[i - 1].width ~/ 2);
      int tRes = Math.max(1, pyramid[i - 1].height ~/ 2);
      pyramid[i] = new SpectrumImage(sRes, tRes, img.samplesPerPixel);

      // Filter four texels from finer level of pyramid
      for (int t = 0, p = 0; t < tRes; ++t) {
        for (int s = 0; s < sRes; ++s, ++p) {
          pyramid[i][p] = (texel(i - 1, 2 * s, 2 * t) +
                           texel(i - 1, 2 * s + 1, 2 * t) +
                           texel(i - 1, 2 * s, 2 * t + 1) +
                           texel(i - 1, 2 * s + 1, 2 * t + 1)) * 0.25;
        }
      }
    }

    // Initialize EWA filter weights if needed
    if (weightLut == null) {
      weightLut = new Float32List(WEIGHT_LUT_SIZE);
      for (int i = 0; i < WEIGHT_LUT_SIZE; ++i) {
        double alpha = 2.0;
        double r2 = i / (WEIGHT_LUT_SIZE - 1);
        weightLut[i] = Math.exp(-alpha * r2) - Math.exp(-alpha);
      }
    }
    if (filename.isNotEmpty && levels > 1) {
      LogInfo('Finished generating MIPMap for $filename');
    }
  }

  texel(int level, int s, int t) {
    assert(level < levels);
    SpectrumImage l = pyramid[level];

    // Compute texel $(s,t)$ accounting for boundary conditions
    switch (wrapMode) {
      case TEXTURE_REPEAT:
        s = s % l.width;
        t = t % l.height;
        break;
      case TEXTURE_CLAMP:
        s = s.clamp(0, l.width - 1);
        t = t.clamp(0, l.height - 1);
        break;
      case TEXTURE_BLACK:
        if (s < 0 || s >= l.width ||
            t < 0 || t >= l.height) {
          return new Spectrum(0.0);
        }
        break;
    }

    return l[t * l.width + s];
  }

  lookup(double s, double t, [double width = 0.0]) {
    // Compute MIPMap level for trilinear filtering
    double level = levels - 1 + Log2(Math.max(width, 1.0e-8));

    // Perform trilinear interpolation at appropriate MIPMap level
    if (level < 0) {
      return triangle(0, s, t);
    } else if (level >= levels - 1) {
      return texel(levels - 1, 0, 0);
    } else {
      int iLevel = level.floor();
      double delta = level - iLevel;
      return triangle(iLevel, s, t) * (1.0 - delta) +
             triangle(iLevel + 1, s, t) * delta;
    }
  }

  lookup2(double s, double t, double ds0, double dt0,
           double ds1, double dt1) {
    if (doTrilinear) {
      Spectrum val = lookup(s, t,
                     2.0 * Math.max(Math.max(ds0.abs(), dt0.abs()),
                                    Math.max(ds1.abs(), dt1.abs())));
      return val;
    }

    // Compute ellipse minor and major axes
    if (ds0 * ds0 + dt0 * dt0 < ds1 * ds1 + dt1 * dt1) {
      double tt = ds0;
      ds0 = ds1;
      ds1 = tt;

      tt = dt0;
      dt0 = dt1;
      dt1 = tt;
    }

    double majorLength = Math.sqrt(ds0 * ds0 + dt0 * dt0);
    double minorLength = Math.sqrt(ds1 * ds1 + dt1 * dt1);

    // Clamp ellipse eccentricity if too large
    if (minorLength * maxAnisotropy < majorLength && minorLength > 0.0) {
      double scale = majorLength / (minorLength * maxAnisotropy);
      ds1 *= scale;
      dt1 *= scale;
      minorLength *= scale;
    }

    if (minorLength == 0.0) {
      var val = triangle(0, s, t);
      return val;
    }

    // Choose level of detail for EWA lookup and perform EWA filtering
    double lod = Math.max(0.0, levels - 1.0 + Log2(minorLength));
    int ilod = lod.floor();

    double d = lod - ilod;
    var val = EWA(ilod, s, t, ds0, dt0, ds1, dt1) * (1.0 - d) +
              EWA(ilod + 1, s, t, ds0, dt0, ds1, dt1) * d;

    return val;
  }

  EWA(int level, double s, double t, double ds0, double dt0,
      double ds1, double dt1) {
    if (level >= levels) {
      return texel(levels - 1, 0, 0);
    }

    // Convert EWA coordinates to appropriate scale for level
    s = s * pyramid[level].width - 0.5;
    t = t * pyramid[level].height - 0.5;
    ds0 *= pyramid[level].width;
    dt0 *= pyramid[level].height;
    ds1 *= pyramid[level].width;
    dt1 *= pyramid[level].height;

    // Compute ellipse coefficients to bound EWA filter region
    double A = dt0 * dt0 + dt1 * dt1 + 1;
    double B = -2.0 * (ds0 * dt0 + ds1 * dt1);
    double C = ds0 * ds0 + ds1 * ds1 + 1;
    double invF = 1.0 / (A * C - B * B * 0.25);
    A *= invF;
    B *= invF;
    C *= invF;

    // Compute the ellipse's $(s,t)$ bounding box in texture space
    double det = -B * B + 4.0 * A * C;
    double invDet = 1.0 / det;
    double uSqrt = Math.sqrt(det * C);
    double vSqrt = Math.sqrt(A * det);
    int s0 = (s - 2.0 * invDet * uSqrt).ceil();
    int s1 = (s + 2.0 * invDet * uSqrt).floor();
    int t0 = (t - 2.0 * invDet * vSqrt).ceil();
    int t1 = (t + 2.0 * invDet * vSqrt).floor();

    // Scan over ellipse bound and compute quadratic equation
    var sum = pyramid[level].samplesPerPixel == 1 ? 0.0 :
              new Spectrum(0.0);

    double sumWts = 0.0;

    for (int it = t0; it <= t1; ++it) {
      double tt = it - t;
      for (int si = s0; si <= s1; ++si) {
        double ss = si - s;
        // Compute squared radius and filter texel if inside ellipse
        double r2 = A * ss * ss + B * ss * tt + C * tt * tt;
        if (r2 < 1.0) {
          double weight = weightLut[Math.min((r2 * WEIGHT_LUT_SIZE),
                                             WEIGHT_LUT_SIZE - 1).toInt()];
          sum += texel(level, si, it) * weight;
          sumWts += weight;
        }
      }
    }

    return sum / sumWts;
  }

  triangle(int level, double s, double t) {
    level = level.clamp(0, levels - 1);
    s = s * pyramid[level].width - 0.5;
    t = t * pyramid[level].height - 0.5;
    int s0 = s.floor();
    int t0 = t.floor();
    double ds = s - s0;
    double dt = t - t0;

    return texel(level, s0, t0) * ((1.0 - ds) * (1.0 - dt)) +
           texel(level, s0, t0 + 1) * ((1.0 - ds) * dt) +
           texel(level, s0 + 1, t0) * (ds * (1.0 - dt)) +
           texel(level, s0 + 1, t0 + 1) * (ds * dt);
  }

  _clamp(f) => f.clamp(0.0, INFINITY);

  List<_ResampleWeight> _resampleWeights(int oldres, int newres) {
    assert(newres >= oldres);
    List<_ResampleWeight> wt = new List<_ResampleWeight>(newres);
    double filterwidth = 2.0;
    for (int i = 0; i < newres; ++i) {
      wt[i] = new _ResampleWeight();

      // Compute image resampling weights for _i_th texel
      double center = (i + 0.5) * oldres / newres;
      wt[i].firstTexel = ((center - filterwidth) + 0.5).floor();
      for (int j = 0; j < 4; ++j) {
        double pos = wt[i].firstTexel + j + 0.5;
        wt[i].weight[j] = Lanczos((pos - center) / filterwidth);
      }

      // Normalize filter weights for texel resampling
      double invSumWts = 1.0 / (wt[i].weight[0] + wt[i].weight[1] +
                                wt[i].weight[2] + wt[i].weight[3]);
      for (int j = 0; j < 4; ++j) {
        wt[i].weight[j] *= invSumWts;
      }
    }

    return wt;
  }

  bool doTrilinear;
  double maxAnisotropy;
  int wrapMode;
  List<SpectrumImage> pyramid;
  int width;
  int height;
  int levels;

  static const int WEIGHT_LUT_SIZE = 128;
  static Float32List weightLut;
}

class _ResampleWeight {
  int firstTexel;
  List<double> weight = [0.0, 0.0, 0.0, 0.0];
}
