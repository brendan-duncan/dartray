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
part of renderers;

class SamplerRenderer extends Renderer {
  SamplerRenderer(this.sampler, this.camera, this.surfaceIntegrator,
                  this.volumeIntegrator,
                  [this.taskNum = 0, this.taskCount = 1]);

  OutputImage render(Scene scene) {
    // Allow integrators to do preprocessing for the scene
    surfaceIntegrator.preprocess(scene, camera, this);

    if (volumeIntegrator != null) {
      volumeIntegrator.preprocess(scene, camera, this);
    }

    // Allocate and initialize sample
    Sampler mainSampler = this.sampler;
    Sample sample = new Sample(mainSampler, surfaceIntegrator,
                               volumeIntegrator, scene);

    Sampler sampler = mainSampler.getSubSampler(taskNum, taskCount);

    RNG rng = new RNG(taskNum);

    // Allocate space for samples and intersections
    int maxSamples = sampler.maximumSampleCount();
    List<Sample> samples = sample.duplicate(maxSamples);

    List<RayDifferential> rays = new List<RayDifferential>(maxSamples);
    List<Spectrum> Ls = new List<Spectrum>(maxSamples);
    List<Spectrum> Ts = new List<Spectrum>(maxSamples);
    List<Intersection> isects = new List<Intersection>(maxSamples);
    for (int i = 0; i < maxSamples; ++i) {
      rays[i] = new RayDifferential();
      Ls[i] = new Spectrum(0.0);
      Ts[i] = new Spectrum(0.0);
      isects[i] = new Intersection();
    }

    // Get samples from [Sampler] and update image
    int sampleCount;
    while ((sampleCount = sampler.getMoreSamples(samples, rng)) > 0) {
      // Generate camera rays and compute radiance along rays
      for (int i = 0; i < sampleCount; ++i) {
        // Find camera ray for sample[i]
        double rayWeight = camera.generateRayDifferential(samples[i], rays[i]);
        rays[i].scaleDifferentials(1.0 / Math.sqrt(sampler.samplesPerPixel));

        // Evaluate radiance along camera ray
        if (rayWeight > 0.0) {
          Ls[i] = Li(scene, rays[i], samples[i], rng, isects[i], Ts[i]) *
                  rayWeight;
        } else {
          Ls[i] = new Spectrum(0.0);
          Ts[i] = new Spectrum(1.0);
        }

        // Issue warning if unexpected radiance value returned
        if (Ls[i].hasNaNs()) {
          LogWarning('Not-a-number radiance value returned '
                'for image sample.  Setting to black.');
          Ls[i] = new Spectrum(0.0);
        } else if (Ls[i].y < -1e-5) {
          LogWarning('Negative luminance value, ${Ls[i].y}, returned'
                'for image sample.  Setting to black.');
          Ls[i] = new Spectrum(0.0);
        } else if (Ls[i].y.isInfinite) {
          LogWarning('Infinite luminance value returned'
                'for image sample.  Setting to black.');
          Ls[i] = new Spectrum(0.0);
        }
      }

      // Report sample results to [Sampler], add contributions to image
      if (sampler.reportResults(samples, rays, Ls, isects, sampleCount)) {
        for (int i = 0; i < sampleCount; ++i) {
          camera.film.addSample(samples[i], Ls[i]);
        }
      }
    }

    // Clean up after [SamplerRenderer] is done with its image region
    camera.film.updateDisplay(sampler.xPixelStart, sampler.yPixelStart,
                              sampler.xPixelEnd + 1, sampler.yPixelEnd + 1);

    return camera.film.writeImage();
  }

  Spectrum Li(Scene scene, RayDifferential ray,
              Sample sample, RNG rng,
              [Intersection isect = null, Spectrum T = null]) {
    if (T == null) {
      T = new Spectrum(0.0);
    }

    if (volumeIntegrator == null) {
      T.set(1.0);
    }

    if (isect == null) {
      isect = new Intersection();
    }

    Spectrum Li = new Spectrum(0.0);

    if (scene.intersect(ray, isect)) {
      Li = surfaceIntegrator.Li(scene, this, ray, isect, sample, rng);
    } else {
      // Handle ray that doesn't intersect any geometry
      for (int i = 0; i < scene.lights.length; ++i) {
        Li += scene.lights[i].Le(ray);
      }
    }

    Spectrum Lvi = volumeIntegrator == null ? new Spectrum(0.0) :
                   volumeIntegrator.Li(scene, this, ray, sample, rng, T);

    return T * Li + Lvi;
  }

  Spectrum transmittance(Scene scene, RayDifferential ray,
                            Sample sample, RNG rng) {
    return volumeIntegrator.transmittance(scene, this, ray, sample, rng);
  }

  int taskNum;
  int taskCount;
  Sampler sampler;
  Camera camera;
  SurfaceIntegrator surfaceIntegrator;
  VolumeIntegrator volumeIntegrator;
}
