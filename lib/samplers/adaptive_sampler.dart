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
part of samplers;

class AdaptiveSampler extends Sampler {
  static const int ADAPTIVE_COMPARE_SHAPE_ID = 0;
  static const int ADAPTIVE_CONTRAST_THRESHOLD = 1;

  AdaptiveSampler(int x, int y, int width, int height, int mins, int maxs,
                  int method, double sopen, double sclose, this.pixels) :
    super(x, y, width, height, sopen, sclose,
          RoundUpPow2(Math.max(mins, maxs))) {
    if (pixels == null) {
      LogSevere('A PixelSampler is required by AdaptiveSampler');
    }
    pixels.setup(x, y, width, height);
    pixelIndex = 0;
    supersamplePixel = false;

    pixels.getPixel(pixelIndex, pixel);

    if (mins > maxs) {
      int t = mins;
      mins = maxs;
      maxs = t;
    }

    if (!IsPowerOf2(mins)) {
      LogWarning('Minimum pixel samples being rounded up to power of 2');
      minSamples = RoundUpPow2(mins);
    } else {
      minSamples = mins;
    }

    if (!IsPowerOf2(maxs)) {
      LogWarning('Maximum pixel samples being rounded up to power of 2');
      maxSamples = RoundUpPow2(maxs);
    } else {
      maxSamples = maxs;
    }

    if (minSamples < 2) {
      LogWarning('Adaptive sampler needs at least two initial pixel samples. '
                 'Using two.');
      minSamples = 2;
    }

    if (minSamples == maxSamples) {
      maxSamples *= 2;
      LogWarning('Adaptive sampler must have more maximum samples than '
                 'minimum. Using $minSamples - $maxSamples');
    }

    pass = 0;
    if (RenderOverrides.SamplingMode() == Sampler.TWO_PASS_SAMPLING ||
        RenderOverrides.SamplingMode() == Sampler.ITERATIVE_SAMPLING) {
      randomSampler = new RandomSampler(x, y, width, height,
                                        sopen, sclose, pixels, 1);
    }
  }

  int roundSize(int size) {
    return RoundUpPow2(size);
  }

  int maximumSampleCount() {
    return maxSamples;
  }

  int getMoreSamples(List<Sample> samples, RNG rng) {
    if (pass == 0 && randomSampler != null) {
      int count = randomSampler.getMoreSamples(samples, rng);
      if (count != 0) {
        return count;
      }
      pass++;
    }

    if (sampleBuf == null) {
      sampleBuf = new Float32List(LDPixelSampleFloatsNeeded(samples[0],
                                                            maxSamples));
    }

    if (supersamplePixel) {
      LDPixelSample(pixel[0], pixel[1], shutterOpen, shutterClose, maxSamples,
                    samples, sampleBuf, rng);
      return maxSamples;
    } else {
      if (pixelIndex >= pixels.numPixels()) {
        return 0;
      }

      LDPixelSample(pixel[0], pixel[1], shutterOpen, shutterClose, minSamples,
                    samples, sampleBuf, rng);

      return minSamples;
    }
  }

  bool reportResults(List<Sample> samples, List<RayDifferential> rays,
                     List<Spectrum> Ls, List<Intersection> isects, int count) {
    if (pass == 0 && randomSampler != null) {
      return true;
    }

    if (supersamplePixel) {
      supersamplePixel = false;
      // Advance to next pixel for sampling
      if (pixelIndex < pixels.numPixels()) {
        pixels.getPixel(++pixelIndex, pixel);
      }
      return true;
    } else if (needsSupersampling(samples, rays, Ls, isects, count)) {
      Stats.SUPERSAMPLE_PIXEL_YES(pixel[0], pixel[1]);
      supersamplePixel = true;
      return false;
    } else {
      Stats.SUPERSAMPLE_PIXEL_NO(pixel[0], pixel[1]);
      // Advance to next pixel for sampling
      if (pixelIndex < pixels.numPixels()) {
        pixels.getPixel(++pixelIndex, pixel);
      }
      return true;
    }
  }

  bool needsSupersampling(List<Sample> samples, List<RayDifferential> rays,
      List<Spectrum> Ls, List<Intersection> isects, int count) {
    switch (method) {
      case ADAPTIVE_COMPARE_SHAPE_ID:
        // See if any shape ids differ within samples
        for (int i = 0; i < count - 1; ++i) {
          if (isects[i].shapeId != isects[i + 1].shapeId ||
              isects[i].primitiveId != isects[i + 1].primitiveId) {
            return true;
          }
        }
        return false;
      case ADAPTIVE_CONTRAST_THRESHOLD:
        // Compare contrast of sample differences to threshold
        double Lavg = 0.0;
        for (int i = 0; i < count; ++i) {
          Lavg += Ls[i].luminance();
        }
        Lavg /= count;
        const double maxContrast = 0.5;
        for (int i = 0; i < count; ++i) {
          if ((Ls[i].luminance() - Lavg).abs() / Lavg > maxContrast) {
            return true;
          }
        }
        return false;
    }

    return false;
  }

  static AdaptiveSampler Create(ParamSet params, int x, int y,
                                int width, int height, Camera camera,
                                PixelSampler pixels) {
    // Initialize common sampler parameters
    int minsamp = params.findOneInt('minsamples', 4);
    int maxsamp = params.findOneInt('maxsamples', 32);

    String m = params.findOneString('method', 'contrast');
    int method = (m == 'contrast') ? ADAPTIVE_CONTRAST_THRESHOLD :
                 (m == 'shapeid') ? ADAPTIVE_COMPARE_SHAPE_ID :
                 -1;

    if (method == -1) {
      LogWarning('Adaptive sampling metric \'$m\' unknown. Using \'contrast\'.');
      method = ADAPTIVE_CONTRAST_THRESHOLD;
    }

    return new AdaptiveSampler(x, y, width, height, minsamp, maxsamp,
                               method, camera.shutterOpen, camera.shutterClose,
                               pixels);
  }

  PixelSampler pixels;
  Int32List pixel = new Int32List(2);
  int pixelIndex;
  int minSamples;
  int maxSamples;
  Float32List sampleBuf;
  int method;
  bool supersamplePixel;
  Sampler randomSampler;
  int pass;
}
