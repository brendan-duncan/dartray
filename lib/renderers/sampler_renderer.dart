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
part of renderers;

/**
 * This is the primary DartRay renderer, able to use any [Sampler],
 * [SurfaceIntegrator] and [VolumeIntegrator].
 */
class SamplerRenderer extends Renderer {
  SamplerRenderer(this.sampler, this.camera, this.surfaceIntegrator,
                  this.volumeIntegrator, [this.taskNum = 0,
                  this.taskCount = 1]);

  /**
   * Render the [scene] from the viewpoint of the [camera].
   */
  Future<OutputImage> render(Scene scene) {
    LogInfo('Starting SamplerRenderer: '
            '${camera.film.xResolution}x${camera.film.yResolution}');

    // Allow integrators to do preprocessing for the scene
    Stats.STARTED_PREPROCESSING();
    surfaceIntegrator.preprocess(scene, camera, this);
    volumeIntegrator.preprocess(scene, camera, this);
    Stats.FINISHED_PREPROCESSING();

    Stats.STARTED_RENDERING();

    // Allocate and initialize sample
    Sample sample = new Sample(sampler, surfaceIntegrator,
                               volumeIntegrator, scene);

    _SamplerRendererTask task = new _SamplerRendererTask(scene, this, camera,
                                                         sampler, sample,
                                                         taskNum, taskCount);

    Completer<OutputImage> completer = new Completer<OutputImage>();

    task.run().then((_) {
      OutputImage out = camera.film.writeImage();
      Stats.FINISHED_RENDERING();
      completer.complete(out);
    });

    return completer.future;
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

class _SamplerRendererTask {
  _SamplerRendererTask(this.scene, this.renderer, this.camera,
                       this.mainSampler, this.origSample,
                       this.taskNum, this.taskCount);

  Future run() {
    Completer completer = new Completer();

    Stats.STARTED_RENDERTASK(taskNum);
    // Get sub-_Sampler_ for _SamplerRendererTask_
    Sampler sampler = mainSampler;
    if (sampler == null) {
      Stats.FINISHED_RENDERTASK(taskNum);
      completer.complete();
      return completer.future;
    }

    Completer renderCompleter = new Completer();

    LogInfo('SamplerRender $taskNum / $taskCount: '
            '[${sampler.left} ${sampler.top} '
            '${sampler.width} ${sampler.height}]');

    // Declare local variables used for rendering loop
    RNG rng = new RNG(taskNum);

    // Allocate space for samples and intersections
    int maxSamples = sampler.maximumSampleCount();
    List<Sample> samples = origSample.duplicate(maxSamples);
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

    while (true) {
      int sampleCount = sampler.getMoreSamples(samples, rng);
      if (sampleCount == 0) {
        renderCompleter.complete();
        break;
      }
      // Generate camera rays and compute radiance along rays
      for (int i = 0; i < sampleCount; ++i) {
        // Find camera ray for _sample[i]_
        Stats.STARTED_GENERATING_CAMERA_RAY(samples[i]);
        double rayWeight = camera.generateRayDifferential(samples[i], rays[i]);
        rays[i].scaleDifferentials(1.0 / Math.sqrt(sampler.samplesPerPixel));
        Stats.FINISHED_GENERATING_CAMERA_RAY(samples[i], rays[i], rayWeight);

        // Evaluate radiance along camera ray
        Stats.STARTED_CAMERA_RAY_INTEGRATION(rays[i], samples[i]);

        if (rayWeight > 0.0) {
          Ls[i] = renderer.Li(scene, rays[i], samples[i], rng, isects[i],
                              Ts[i]) * rayWeight;
        } else {
          Ls[i] = new Spectrum(0.0);
          Ts[i] = new Spectrum(1.0);
        }

        // Issue warning if unexpected radiance value returned
        if (Ls[i].hasNaNs()) {
          LogWarning('Not-a-number radiance value returned '
                     'for image sample. Setting to black.');
          Ls[i] = new Spectrum(0.0);
        } else if (Ls[i].luminance() < -1e-5) {
          LogWarning('Negative luminance value, ${Ls[i].luminance()}, returned'
                     'for image sample. Setting to black.');
          Ls[i] = new Spectrum(0.0);
        } else if (Ls[i].luminance().isInfinite) {
          LogWarning('Infinite luminance value returned'
                     'for image sample. Setting to black.');
          Ls[i] = new Spectrum(0.0);
        }

        Stats.FINISHED_CAMERA_RAY_INTEGRATION(rays[i], samples[i], Ls[i]);
      }

      // Report sample results to [Sampler], add contributions to image
      if (sampler.reportResults(samples, rays, Ls, isects, sampleCount)) {
        for (int i = 0; i < sampleCount; ++i) {
          Stats.STARTED_ADDING_IMAGE_SAMPLE(samples[i], rays[i], Ls[i], Ts[i]);
          camera.film.addSample(samples[i], Ls[i]);
          Stats.FINISHED_ADDING_IMAGE_SAMPLE();
        }
      }
    }

    renderCompleter.future.then((_) {
      // Clean up after [SamplerRenderer] is done with its image region
      camera.film.updateDisplay(sampler.left, sampler.top,
                                sampler.width, sampler.height);

      Stats.FINISHED_RENDERTASK(taskNum);
      completer.complete();
    });

    return completer.future;
  }

  Scene scene;
  Renderer renderer;
  Camera camera;
  Sampler mainSampler;
  Sample origSample;
  int taskNum;
  int taskCount;
}
