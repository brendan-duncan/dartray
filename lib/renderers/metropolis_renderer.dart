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

class MetropolisRenderer extends Renderer {
  MetropolisRenderer(this.nPixelSamples, this.nBootstrap,
                     this.nDirectPixelSamples, double lsp,
                     bool doDirectSeparately, this.maxConsecutiveRejects,
                     this.maxDepth, this.camera, this.bidirectional,
                     [this.taskNum = 0, this.taskCount = 1]) {
    largeStepsPerPixel = Math.max(1, RoundUpPow2((lsp * nPixelSamples).toInt()));
    if (largeStepsPerPixel >= nPixelSamples && largeStepsPerPixel > 1) {
      largeStepsPerPixel ~/= 2;
    }

    assert(largeStepsPerPixel >= 1 && largeStepsPerPixel < nPixelSamples);

    if ((nPixelSamples % largeStepsPerPixel) != 0) {
      int origPixelSamples = nPixelSamples;
      nPixelSamples += largeStepsPerPixel - (nPixelSamples % largeStepsPerPixel);
      LogWarning('Rounding up to $nPixelSamples Metropolis samples '
                 'per pixel (from $origPixelSamples)');
    }

    directLighting = doDirectSeparately ?
        Plugin.surfaceIntegrator('directlighting')(new ParamSet.fromJson({
          'int maxdepth': [maxDepth],
          'string strategy': ['all'] })) :
        null;
  }

  Future<OutputImage> render(Scene scene) {
    LogInfo('Starting MetropolisRenderer: '
                '${camera.film.xResolution}x${camera.film.yResolution}');
    Stats.MLT_STARTED_RENDERING();

    if (scene.lights.length > 0) {
      List<int> extent = [0, 0, 0, 0];
      camera.film.getPixelExtent(extent);

      double t0 = camera.shutterOpen;
      double t1 = camera.shutterClose;
      Distribution1D lightDistribution =
          Integrator.ComputeLightSamplingCDF(scene);

      if (directLighting != null) {
        Stopwatch t = new Stopwatch()..start();
        LogInfo('Metropolis: Starting Direct Lighting Render');
        Stats.MLT_STARTED_DIRECTLIGHTING();
        // Compute direct lighting before Metropolis light transport
        if (nDirectPixelSamples > 0) {
          LowDiscrepancySampler sampler =
              new LowDiscrepancySampler(extent[0], extent[2],
                                        extent[1] - extent[0],
                                        extent[3] - extent[2],
                                        t0, t1, new TilePixelSampler(),
                                        nDirectPixelSamples);
          Sample sample = new Sample(sampler, directLighting, null, scene);

          var task = new _SamplerRendererTask(scene, this, camera,
                                              sampler, sample,
                                              taskNum, taskCount);
          task.run();
        }

        camera.film.writeImage();
        Stats.MLT_FINISHED_DIRECTLIGHTING();
        LogInfo('Metropolis:Finished Direct Lighting Render: ${t.elapsed}');
      }

      Stopwatch t = new Stopwatch()..start();
      LogInfo('Metropolis: Starting Bootstrap: $nBootstrap');

      // Take initial set of samples to compute $b$
      Stats.MLT_STARTED_BOOTSTRAPPING(nBootstrap);
      RNG rng = new RNG(taskNum);

      GetSubWindow((extent[1] - extent[0]) + 1, (extent[3] - extent[2]) + 1,
                               taskNum, taskCount, extent);
      LogInfo('RENDER EXTENT $extent');

      List<_PathVertex> cameraPath = new List<_PathVertex>(maxDepth);
      List<_PathVertex> lightPath = new List<_PathVertex>(maxDepth);
      for (int i = 0; i < maxDepth; ++i) {
        cameraPath[i] = new _PathVertex();
        lightPath[i] = new _PathVertex();
      }

      double sumI = 0.0;
      Float32List bootstrapI = new Float32List(nBootstrap);

      _MLTSample sample = new _MLTSample(maxDepth);

      for (int i = 0; i < nBootstrap; ++i) {
        // Generate random sample and path radiance for MLT bootstrapping
        double x = Lerp(rng.randomFloat(), extent[0], extent[1]);
        double y = Lerp(rng.randomFloat(), extent[2], extent[3]);

        LargeStep(rng, sample, maxDepth, x, y, t0, t1, bidirectional);

        Spectrum L = pathL(sample, scene, camera, lightDistribution,
                           cameraPath, lightPath, rng);

        // Compute contribution for random sample for MLT bootstrapping
        double I = L.luminance();
        sumI += I;
        bootstrapI[i] = I;
      }

      double b = sumI / nBootstrap;

      LogInfo('Metropolis: Finished Bootstrap: ${t.elapsed}');
      Stats.MLT_FINISHED_BOOTSTRAPPING(b);

      LogInfo("MLT computed b = $b");

      // Select initial sample from bootstrap samples
      double contribOffset = rng.randomFloat() * sumI;
      rng.seed(0);
      sumI = 0.0;
      _MLTSample initialSample = new _MLTSample(maxDepth);

      for (int i = 0; i < nBootstrap; ++i) {
        double x = Lerp(rng.randomFloat(), extent[0], extent[1]);
        double y = Lerp(rng.randomFloat(), extent[2], extent[3]);

        LargeStep(rng, initialSample, maxDepth, x, y, t0, t1, bidirectional);

        sumI += bootstrapI[i];

        if (sumI > contribOffset) {
          break;
        }
      }

      // Launch tasks to generate Metropolis samples
      int nTasks = largeStepsPerPixel;
      int largeStepRate = nPixelSamples ~/ largeStepsPerPixel;
      LogInfo("MLT running $nTasks tasks, large step rate $largeStepRate");

      List<int> scramble = [ rng.randomUint(), rng.randomUint() ];

      for (int i = 0; i < nTasks; ++i) {
        Stopwatch t = new Stopwatch()..start();
        LogInfo('Metropolis: Task ${i + 1} / $nTasks');
        List<double> d = [0.0, 0.0];
        Sample02(0, scramble, d);

        _MLTTask task = new _MLTTask(i, d[0], d[1], extent[0], extent[1],
                                     extent[2], extent[3], t0, t1, b,
                                     initialSample, scene, camera, this,
                                     lightDistribution);
        task.run();
        LogInfo('Metropolis: Finished Task ${i + 1} / $nTasks: ${t.elapsed}');
      }
    }

    OutputImage out = camera.film.writeImage();

    Stats.MLT_FINISHED_RENDERING();
    Completer<OutputImage> c = new Completer<OutputImage>();
    c.complete(out);

    return c.future;
  }

  Spectrum Li(Scene scene, RayDifferential ray,
              Sample sample, RNG rng, [Intersection isect, Spectrum T]) {
    Intersection localIsect = new Intersection();
    if (isect == null) {
      isect = localIsect;
    }

    Spectrum Lo = new Spectrum(0.0);

    if (scene.intersect(ray, isect)) {
      Lo = directLighting.Li(scene, this, ray, isect, sample, rng);
    } else {
      // Handle ray that doesn't intersect any geometry
      for (int i = 0; i < scene.lights.length; ++i) {
        Lo += scene.lights[i].Le(ray);
      }
    }

    return Lo;
  }

  Spectrum transmittance(Scene scene, RayDifferential ray, Sample sample,
                         RNG rng) {
    return new Spectrum(1.0);
  }

  Spectrum pathL(_MLTSample sample, Scene scene, Camera camera,
                 Distribution1D lightDistribution,
                 List<_PathVertex> cameraPath, List<_PathVertex> lightPath,
                 RNG rng) {
    // Generate camera path from camera path samples
    Stats.STARTED_GENERATING_CAMERA_RAY(sample.cameraSample);

    RayDifferential cameraRay = new RayDifferential();
    double cameraWt = camera.generateRayDifferential(sample.cameraSample,
                                                     cameraRay);
    cameraRay.scaleDifferentials(1.0 / Math.sqrt(nPixelSamples));

    Stats.FINISHED_GENERATING_CAMERA_RAY(sample.cameraSample, cameraRay,
                                         cameraWt);

    RayDifferential escapedRay = new RayDifferential();
    Spectrum escapedAlpha = new Spectrum();

    int cameraLength = GeneratePath(cameraRay, new Spectrum(cameraWt), scene,
                                    sample.cameraPathSamples, cameraPath,
                                    escapedRay, escapedAlpha);

    if (!bidirectional) {
      // Compute radiance along path using path tracing
      return Lpath(scene, cameraPath, cameraLength, sample.lightingSamples,
                   rng, sample.cameraSample.time, lightDistribution,
                   escapedRay, escapedAlpha);
    } else {
      // Sample light ray and apply bidirectional path tracing

      // Choose light and sample ray to start light path
      Stats.MLT_STARTED_SAMPLE_LIGHT_FOR_BIDIR();

      List<double> lightPdf = [0.0];
      List<double> lightRayPdf = [0.0];

      int lightNum = lightDistribution.sampleDiscrete(sample.lightNumSample,
                                                      lightPdf);

      Light light = scene.lights[lightNum];
      Ray lightRay = new Ray();
      Normal Nl = new Normal();
      LightSample lrs = new LightSample(sample.lightRaySamples[0],
                                        sample.lightRaySamples[1],
                                        sample.lightRaySamples[2]);

      Spectrum lightWt = light.sampleL(scene, lrs, sample.lightRaySamples[3],
                                       sample.lightRaySamples[4],
                                       sample.cameraSample.time, lightRay,
                                       Nl, lightRayPdf);

      Stats.MLT_FINISHED_SAMPLE_LIGHT_FOR_BIDIR();

      if (lightWt.isBlack() || lightRayPdf[0] == 0.0) {
        // Compute radiance along path using path tracing
        return Lpath(scene, cameraPath, cameraLength,
                     sample.lightingSamples, rng, sample.cameraSample.time,
                     lightDistribution, escapedRay, escapedAlpha);
      } else {
        // Compute radiance along paths using bidirectional path tracing
        lightWt *= Vector.AbsDot(Normal.Normalize(Nl), lightRay.direction) /
                   (lightPdf[0] * lightRayPdf[0]);

        int lightLength = GeneratePath(new RayDifferential.fromRay(lightRay),
                                       lightWt, scene, sample.lightPathSamples,
                                       lightPath, null, null);

        return Lbidir(scene, cameraPath, cameraLength, lightPath, lightLength,
                      sample.lightingSamples, rng, sample.cameraSample.time,
                      lightDistribution, escapedRay, escapedAlpha);
      }
    }
  }

  Spectrum Lpath(Scene scene, List<_PathVertex> cameraPath,
                 int cameraPathLength, List<_LightingSample> samples,
                 RNG rng, double time, Distribution1D lightDistribution,
                 RayDifferential escapedRay, Spectrum escapedAlpha) {
    Stats.MLT_STARTED_LPATH();

    Spectrum L = new Spectrum(0.0);
    bool previousSpecular = true;
    bool allSpecular = true;

    for (int i = 0; i < cameraPathLength; ++i) {
      // Initialize basic variables for camera path vertex
      _PathVertex vc = cameraPath[i];
      Point pc = vc.bsdf.dgShading.p;
      Normal nc = vc.bsdf.dgShading.nn;

      // Add emitted light from vertex if appropriate
      if (previousSpecular && (directLighting == null || !allSpecular)) {
        L += vc.alpha * vc.isect.Le(vc.wPrev);
      }

      // Compute direct illumination for Metropolis path vertex
      Spectrum Ld = new Spectrum(0.0);

      if (directLighting == null || !allSpecular) {
        // Choose light and call _EstimateDirect()_ for Metropolis vertex
        _LightingSample ls = samples[i];
        List<double> lightPdf = [0.0];
        int lightNum = lightDistribution.sampleDiscrete(ls.lightNum, lightPdf);
        Light light = scene.lights[lightNum];

        Stats.MLT_STARTED_ESTIMATE_DIRECT();

        Ld = vc.alpha *
             Integrator.EstimateDirect(scene, this, light, pc, nc,
                                       vc.wPrev, vc.isect.rayEpsilon,
                                       time, vc.bsdf, rng,
                                       ls.lightSample, ls.bsdfSample,
                                       BSDF_ALL & ~BSDF_SPECULAR) / lightPdf[0];

        Stats.MLT_FINISHED_ESTIMATE_DIRECT();
      }

      previousSpecular = vc.specularBounce;

      allSpecular = allSpecular && previousSpecular;

      L += Ld;
    }

    // Add contribution of escaped ray, if any
    if (!escapedAlpha.isBlack() && previousSpecular &&
        (directLighting == null || !allSpecular)) {
      for (int i = 0; i < scene.lights.length; ++i) {
        L += escapedAlpha * scene.lights[i].Le(escapedRay);
      }
    }

    Stats.MLT_FINISHED_LPATH();

    return L;
  }

  Spectrum Lbidir(Scene scene,
                  List<_PathVertex> cameraPath, int cameraPathLength,
                  List<_PathVertex> lightPath, int lightPathLength,
                  List<_LightingSample> samples, RNG rng, double time,
                  Distribution1D lightDistribution,
                  RayDifferential escapedRay, Spectrum escapedAlpha) {
    Stats.MLT_STARTED_LBIDIR();

    Spectrum L = new Spectrum(0.0);
    bool previousSpecular = true;
    bool allSpecular = true;

    // Compute number of specular vertices for each path length
    int nVerts = cameraPathLength + lightPathLength + 2;
    Uint32List nSpecularVertices = new Uint32List(nVerts);

    for (int i = 0; i < cameraPathLength; ++i) {
      for (int j = 0; j < lightPathLength; ++j) {
        if (cameraPath[i].specularBounce || lightPath[j].specularBounce) {
          ++nSpecularVertices[i + j + 2];
        }
      }
    }

    for (int i = 0; i < cameraPathLength; ++i) {
      // Initialize basic variables for camera path vertex
      _PathVertex vc = cameraPath[i];
      Point pc = vc.bsdf.dgShading.p;
      Normal nc = vc.bsdf.dgShading.nn;

      // Compute reflected light at camera path vertex

      // Add emitted light from vertex if appropriate
      if (previousSpecular && (directLighting == null || !allSpecular)) {
        L += vc.alpha * vc.isect.Le(vc.wPrev);
      }

      // Compute direct illumination for Metropolis path vertex
      Spectrum Ld = new Spectrum(0.0);
      if (directLighting == null || !allSpecular) {
        // Choose light and call _EstimateDirect()_ for Metropolis vertex
        _LightingSample ls = samples[i];
        List<double> lightPdf = [0.0];
        int lightNum = lightDistribution.sampleDiscrete(ls.lightNum, lightPdf);
        Light light = scene.lights[lightNum];

        Stats.MLT_STARTED_ESTIMATE_DIRECT();

        Ld = vc.alpha *
             Integrator.EstimateDirect(scene, this, light, pc, nc, vc.wPrev,
                                       vc.isect.rayEpsilon, time, vc.bsdf, rng,
                                       ls.lightSample, ls.bsdfSample,
                                       BSDF_ALL & ~BSDF_SPECULAR) / lightPdf[0];

        Stats.MLT_FINISHED_ESTIMATE_DIRECT();
      }

      previousSpecular = vc.specularBounce;
      allSpecular = allSpecular && previousSpecular;
      L += Ld / (i + 1 - nSpecularVertices[i + 1]);

      if (!vc.specularBounce) {
        // Loop over light path vertices and connect to camera vertex
        for (int j = 0; j < lightPathLength; ++j) {
          _PathVertex vl = lightPath[j];
          Point pl = vl.bsdf.dgShading.p;
          Normal nl = vl.bsdf.dgShading.nn;
          if (!vl.specularBounce) {
            // Compute contribution between camera and light vertices
            Vector w = Vector.Normalize(pl - pc);
            Spectrum fc = vc.bsdf.f(vc.wPrev, w) * (1 + vc.nSpecularComponents);
            Spectrum fl = vl.bsdf.f(-w, vl.wPrev) * (1 + vl.nSpecularComponents);

            if (fc.isBlack() || fl.isBlack()) {
              continue;
            }

            Ray r = new Ray(pc, pl - pc, 1.0e-3, 0.999, time);

            if (!scene.intersectP(r)) {
              // Compute weight for bidirectional path, _pathWt_
              double pathWt = 1.0 / (i + j + 2 - nSpecularVertices[i + j + 2]);
              double G = Vector.AbsDot(nc, w) * Vector.AbsDot(nl, w) /
                         Vector.DistanceSquared(pl, pc);
              L += (vc.alpha * fc * G * fl * vl.alpha) * pathWt;
            }
          }
        }
      }
    }

    // Add contribution of escaped ray, if any
    if (!escapedAlpha.isBlack() && previousSpecular &&
        (directLighting == null || !allSpecular)) {
      for (int i = 0; i < scene.lights.length; ++i) {
        L += escapedAlpha * scene.lights[i].Le(escapedRay);
      }
    }

    Stats.MLT_FINISHED_LBIDIR();

    return L;
  }

  static void LargeStep(RNG rng, _MLTSample sample, int maxDepth,
                        double x, double y, double t0, double t1,
                        bool bidirectional) {
    // Do large step mutation of _cameraSample_
    sample.cameraSample.imageX = x;
    sample.cameraSample.imageY = y;
    sample.cameraSample.time = Lerp(rng.randomFloat(), t0, t1);
    sample.cameraSample.lensU = rng.randomFloat();
    sample.cameraSample.lensV = rng.randomFloat();

    for (int i = 0; i < maxDepth; ++i) {
      // Apply large step to i'th camera PathSample
      _PathSample cps = sample.cameraPathSamples[i];
      cps.bsdfSample.uComponent = rng.randomFloat();
      cps.bsdfSample.uDir[0] = rng.randomFloat();
      cps.bsdfSample.uDir[1] = rng.randomFloat();
      cps.rrSample = rng.randomFloat();

      // Apply large step to i'th LightingSample
      _LightingSample ls = sample.lightingSamples[i];
      ls.bsdfSample.uComponent = rng.randomFloat();
      ls.bsdfSample.uDir[0] = rng.randomFloat();
      ls.bsdfSample.uDir[1] = rng.randomFloat();
      ls.lightNum = rng.randomFloat();
      ls.lightSample.uComponent = rng.randomFloat();
      ls.lightSample.uPos[0] = rng.randomFloat();
      ls.lightSample.uPos[1] = rng.randomFloat();
    }

    if (bidirectional) {
      // Apply large step to bidirectional light samples
      sample.lightNumSample = rng.randomFloat();

      for (int i = 0; i < 5; ++i) {
        sample.lightRaySamples[i] = rng.randomFloat();
      }

      for (int i = 0; i < maxDepth; ++i) {
        // Apply large step to i'th light PathSample
        _PathSample lps = sample.lightPathSamples[i];
        lps.bsdfSample.uComponent = rng.randomFloat();
        lps.bsdfSample.uDir[0] = rng.randomFloat();
        lps.bsdfSample.uDir[1] = rng.randomFloat();
        lps.rrSample = rng.randomFloat();
      }
    }
  }

  static double Mutate(RNG rng, double v, [num min = 0, num max = 1]) {
    if (min == max) {
      v = min;
      return v;
    }
    assert(min < max);

    double a = 1.0 / 1024.0;
    double b = 1.0 / 64.0;

    final double logRatio = -Math.log(b / a);
    double delta = (max - min) * b * Math.exp(logRatio * rng.randomFloat());
    if (rng.randomFloat() < 0.5) {
      v += delta;
      if (v >= max) {
        v = min + (v - max);
      }
    } else {
      v -= delta;
      if (v < min) {
        v = max - (min - v);
      }
    }

    if (v < min || v >= max) {
      v = min;
    }

    return v;
  }

  static void SmallStep(RNG rng, _MLTSample sample, int maxDepth,
                        int x0, int x1, int y0, int y1, double t0, double t1,
                        bool bidirectional) {
    sample.cameraSample.imageX = Mutate(rng, sample.cameraSample.imageX, x0, x1);
    sample.cameraSample.imageY = Mutate(rng, sample.cameraSample.imageY, y0, y1);
    sample.cameraSample.time = Mutate(rng, sample.cameraSample.time, t0, t1);
    sample.cameraSample.lensU = Mutate(rng, sample.cameraSample.lensU);
    sample.cameraSample.lensV = Mutate(rng, sample.cameraSample.lensV);

    // Apply small step mutation to camera, lighting, and light samples
    for (int i = 0; i < maxDepth; ++i) {
      // Apply small step to $i$th camera _PathSample_
      _PathSample eps = sample.cameraPathSamples[i];
      eps.bsdfSample.uComponent = Mutate(rng, eps.bsdfSample.uComponent);
      eps.bsdfSample.uDir[0] = Mutate(rng, eps.bsdfSample.uDir[0]);
      eps.bsdfSample.uDir[1] = Mutate(rng, eps.bsdfSample.uDir[1]);
      eps.rrSample = Mutate(rng, eps.rrSample);

      // Apply small step to $i$th _LightingSample_
      _LightingSample ls = sample.lightingSamples[i];
      ls.bsdfSample.uComponent = Mutate(rng, ls.bsdfSample.uComponent);
      ls.bsdfSample.uDir[0] = Mutate(rng, ls.bsdfSample.uDir[0]);
      ls.bsdfSample.uDir[1] = Mutate(rng, ls.bsdfSample.uDir[1]);
      ls.lightNum = Mutate(rng, ls.lightNum);
      ls.lightSample.uComponent = Mutate(rng, ls.lightSample.uComponent);
      ls.lightSample.uPos[0] = Mutate(rng, ls.lightSample.uPos[0]);
      ls.lightSample.uPos[1] = Mutate(rng, ls.lightSample.uPos[1]);
    }

    if (bidirectional) {
      sample.lightNumSample = Mutate(rng, sample.lightNumSample);
      for (int i = 0; i < 5; ++i) {
        sample.lightRaySamples[i] = Mutate(rng, sample.lightRaySamples[i]);
      }
      for (int i = 0; i < maxDepth; ++i) {
        // Apply small step to $i$th light _PathSample_
        _PathSample lps = sample.lightPathSamples[i];
        lps.bsdfSample.uComponent = Mutate(rng, lps.bsdfSample.uComponent);
        lps.bsdfSample.uDir[0] = Mutate(rng, lps.bsdfSample.uDir[0]);
        lps.bsdfSample.uDir[1] = Mutate(rng, lps.bsdfSample.uDir[1]);
        lps.rrSample = Mutate(rng, lps.rrSample);
      }
    }
  }

  static int GeneratePath(RayDifferential r, Spectrum a, Scene scene,
                          List<_PathSample> samples, List<_PathVertex> path,
                          RayDifferential escapedRay, Spectrum escapedAlpha) {
    Stats.MLT_STARTED_GENERATE_PATH();
    RayDifferential ray = r;
    Spectrum alpha = a;
    if (escapedAlpha != null) {
      escapedAlpha.set(0.0);
    }

    int length = 0;
    for (; length < samples.length; ++length) {
      // Try to generate next vertex of ray path
      _PathVertex v = path[length];
      if (!scene.intersect(ray, v.isect)) {
        // Handle ray that leaves the scene during path generation
        if (escapedAlpha != null) {
          escapedAlpha.copy(alpha);
        }
        if (escapedRay != null) {
          escapedRay.copy(ray);
        }
        break;
      }

      // Record information for current path vertex
      v.alpha = alpha;
      BSDF bsdf = v.isect.getBSDF(ray);
      v.bsdf = bsdf;
      v.wPrev = -ray.direction;

      // Sample direction for outgoing Metropolis path direction
      List<double> pdf = [0.0];
      List<int> flags = [0];
      Spectrum f = bsdf.sample_f(-ray.direction, v.wNext,
                                 samples[length].bsdfSample,
                                 pdf, BSDF_ALL, flags);
      v.specularBounce = (flags[0] & BSDF_SPECULAR) != 0;
      v.nSpecularComponents = bsdf.numComponents(BSDF_SPECULAR |
                                       BSDF_REFLECTION | BSDF_TRANSMISSION);
      if (f.isBlack() || pdf[0] == 0.0) {
        Stats.MLT_FINISHED_GENERATE_PATH();
        return length + 1;
      }

      // Terminate path with RR or prepare for finding next vertex
      Point p = bsdf.dgShading.p;
      Normal n = bsdf.dgShading.nn;
      Spectrum pathScale = f * Vector.AbsDot(v.wNext, n) / pdf[0];
      double rrSurviveProb = Math.min(1.0, pathScale.luminance());
      if (samples[length].rrSample > rrSurviveProb) {
        Stats.MLT_FINISHED_GENERATE_PATH();
        return length + 1;
      }

      alpha *= pathScale / rrSurviveProb;

      //alpha *= renderer.tansmittance(scene, ray, null, rng);
      ray = new RayDifferential.child(p, v.wNext, ray, v.isect.rayEpsilon);
    }

    Stats.MLT_FINISHED_GENERATE_PATH();
    return length;
  }


  static MetropolisRenderer Create(ParamSet params, Camera camera,
                                   [int taskNum = 0, int taskCount = 1]) {
    double largeStepProb = params.findOneFloat('largestepprobability', 0.25);
    int perPixelSamples = params.findOneInt('samplesperpixel', 100);
    int nBootstrap = params.findOneInt('bootstrapsamples', 100000);
    int nDirectPixelSamples = params.findOneInt('directsamples', 4);
    bool doDirectSeparately = params.findOneBool('dodirectseparately', true);
    int mr = params.findOneInt('maxconsecutiverejects', 512);
    int md = params.findOneInt('maxdepth', 7);
    bool doBidirectional = params.findOneBool('bidirectional', true);

    if (RenderOverrides.QuickRender()) {
      perPixelSamples = Math.max(1, perPixelSamples ~/ 4);
      nBootstrap = Math.max(1, nBootstrap ~/ 4);
      nDirectPixelSamples = Math.max(1, nDirectPixelSamples ~/ 4);
    }

    return new MetropolisRenderer(perPixelSamples, nBootstrap,
                                  nDirectPixelSamples, largeStepProb,
                                  doDirectSeparately, mr, md, camera,
                                  doBidirectional, taskNum, taskCount);
  }

  int taskNum;
  int taskCount;
  Camera camera;
  bool bidirectional;
  int nDirectPixelSamples;
  int nPixelSamples;
  int maxDepth;
  int largeStepsPerPixel;
  int nBootstrap;
  int maxConsecutiveRejects;
  SurfaceIntegrator directLighting;
}

class _PathSample {
  BSDFSample bsdfSample;
  double rrSample;

  _PathSample() :
    bsdfSample = new BSDFSample(),
    rrSample = 0.0;

  _PathSample.from(_PathSample other) :
    bsdfSample = new BSDFSample.from(other.bsdfSample),
    rrSample = other.rrSample;
}

class _LightingSample {
  BSDFSample bsdfSample;
  double lightNum;
  LightSample lightSample;

  _LightingSample() :
    bsdfSample = new BSDFSample(),
    lightNum = 0.0,
    lightSample = new LightSample();

  _LightingSample.from(_LightingSample other) :
    bsdfSample = new BSDFSample.from(other.bsdfSample),
    lightNum = other.lightNum,
    lightSample = new LightSample.from(other.lightSample);
}

class _PathVertex {
  Intersection isect;
  Vector wPrev;
  Vector wNext;
  BSDF bsdf;
  bool specularBounce;
  int nSpecularComponents;
  Spectrum alpha;

  _PathVertex() :
    isect = new Intersection(),
    wPrev = new Vector(),
    wNext = new Vector(),
    alpha = new Spectrum(0.0);
}

class _MLTSample {
  _MLTSample(int maxLength) :
    cameraSample = new CameraSample(),
    lightNumSample = 0.0,
    cameraPathSamples = new List<_PathSample>(maxLength),
    lightPathSamples = new List<_PathSample>(maxLength),
    lightingSamples = new List<_LightingSample>(maxLength) {
    for (int i = 0; i < maxLength; ++i) {
      cameraPathSamples[i] = new _PathSample();
      lightPathSamples[i] = new _PathSample();
      lightingSamples[i] = new _LightingSample();
    }
  }

  _MLTSample.from(_MLTSample other) :
    cameraSample = new CameraSample.from(other.cameraSample),
    lightNumSample = other.lightNumSample,
    lightRaySamples = new Float32List.fromList(other.lightRaySamples),
    cameraPathSamples = new List<_PathSample>.from(other.cameraPathSamples),
    lightPathSamples = new List<_PathSample>.from(other.lightPathSamples),
    lightingSamples = new List<_LightingSample>.from(other.lightingSamples) {
    for (int i = 0; i < cameraPathSamples.length; ++i) {
      cameraPathSamples[i] = new _PathSample.from(cameraPathSamples[i]);
    }
    for (int i = 0; i < lightPathSamples.length; ++i) {
      lightPathSamples[i] = new _PathSample.from(lightPathSamples[i]);
    }
    for (int i = 0; i < lightingSamples.length; ++i) {
      lightingSamples[i] = new _LightingSample.from(lightingSamples[i]);
    }
  }

  CameraSample cameraSample;
  double lightNumSample;
  Float32List lightRaySamples = new Float32List(5);
  List<_PathSample> cameraPathSamples;
  List<_PathSample> lightPathSamples;
  List<_LightingSample> lightingSamples;
}

class _MLTTask {
  _MLTTask(this.taskNum, this.dx, this.dy, this.x0, this.x1, this.y0, this.y1,
           this.t0, this.t1, this.b, this.initialSample,
           this.scene, this.camera, this.renderer,
           this.lightDistribution) {
    currentPixelSample = 0;
  }

  void run() {
    Stats.MLT_STARTED_MLT_TASK(this);
    // Declare basic _MLTTask_ variables and prepare for sampling
    Stats.MLT_STARTED_TASK_INIT();
    int nPixels = (x1 - x0) * (y1 - y0);
    int nPixelSamples = renderer.nPixelSamples;
    int largeStepRate = nPixelSamples ~/ renderer.largeStepsPerPixel;
    assert(largeStepRate > 1);
    int nTaskSamples = nPixels * largeStepRate;
    int consecutiveRejects = 0;

    int ntf = taskNum + 1;
    int totalSamples = nPixels * nPixelSamples;
    double splatScale = totalSamples / (ntf * nTaskSamples);

    camera.film.splatScale = splatScale;

    // Declare variables for storing and computing MLT samples
    RNG rng = new RNG(taskNum);
    List<_PathVertex> cameraPath = new List<_PathVertex>(renderer.maxDepth);
    List<_PathVertex> lightPath = new List<_PathVertex>(renderer.maxDepth);
    for (int i = 0; i < renderer.maxDepth; ++i) {
      cameraPath[i] = new _PathVertex();
      lightPath[i] = new _PathVertex();
    }

    List<_MLTSample> samples = new List<_MLTSample>(2);
    for (int i = 0; i < 2; ++i) {
      samples[i] = new _MLTSample(renderer.maxDepth);
    }

    List<Spectrum> L = new List<Spectrum>(2);
    Float32List I = new Float32List(2);
    int current = 0;
    int proposed = 1;

    // Compute _L[current]_ for initial sample
    samples[current] = new _MLTSample.from(initialSample);
    L[current] = renderer.pathL(initialSample, scene, camera,
                                lightDistribution, cameraPath, lightPath,
                                rng);
    I[current] = L[current].luminance();

    // Compute randomly permuted table of pixel indices for large steps
    int pixelNumOffset = 0;
    Uint32List largeStepPixelNum = new Uint32List(nPixels);

    for (int i = 0; i < nPixels; ++i) {
      largeStepPixelNum[i] = i;
    }

    Shuffle(largeStepPixelNum, 0, nPixels, 1, rng);

    Stats.MLT_FINISHED_TASK_INIT();

    for (int s = 0; s < nTaskSamples; ++s) {
      // Compute proposed mutation to current sample
      Stats.MLT_STARTED_MUTATION();
      samples[proposed] = new _MLTSample.from(samples[current]);
      bool largeStep = ((s % largeStepRate) == 0);

      if (largeStep) {
        int x = x0 + largeStepPixelNum[pixelNumOffset] % (x1 - x0);
        int y = y0 + largeStepPixelNum[pixelNumOffset] ~/ (x1 - x0);

        MetropolisRenderer.LargeStep(rng, samples[proposed], renderer.maxDepth,
                                     x + dx, y + dy, t0, t1,
                                     renderer.bidirectional);

        ++pixelNumOffset;
      } else {
        MetropolisRenderer.SmallStep(rng, samples[proposed], renderer.maxDepth,
                                     x0, x1, y0, y1, t0, t1,
                                     renderer.bidirectional);
      }

      Stats.MLT_FINISHED_MUTATION();

      // Compute contribution of proposed sample
      L[proposed] = renderer.pathL(samples[proposed], scene, camera,
                                   lightDistribution, cameraPath, lightPath,
                                   rng);
      I[proposed] = L[proposed].luminance();

      // Compute acceptance probability for proposed sample
      double a = Math.min(1.0, I[proposed] / I[current]);

      // Splat current and proposed samples to _Film_
      Stats.MLT_STARTED_SAMPLE_SPLAT();

      if (I[current] > 0.0) {
        if ((1.0 / I[current]).isFinite) {
          Spectrum contrib = L[current] * (b / nPixelSamples) / I[current];
          camera.film.splat(samples[current].cameraSample, contrib * (1.0 - a));
        }
      }

      if (I[proposed] > 0.0) {
        if ((1.0 / I[proposed]).isFinite) {
          Spectrum contrib = L[proposed] * (b / nPixelSamples) / I[proposed];
          camera.film.splat(samples[proposed].cameraSample, contrib * a);
        }
      }

      Stats.MLT_FINISHED_SAMPLE_SPLAT();

      // Randomly accept proposed path mutation (or not)
      if (consecutiveRejects >= renderer.maxConsecutiveRejects ||
          rng.randomFloat() < a) {
        Stats.MLT_ACCEPTED_MUTATION(a, samples[current], samples[proposed]);
        current ^= 1;
        proposed ^= 1;
        consecutiveRejects = 0;
      } else {
        Stats.MLT_REJECTED_MUTATION(a, samples[current], samples[proposed]);
        ++consecutiveRejects;
      }
    }

    assert(pixelNumOffset == nPixels);

    // Update display for recently computed Metropolis samples
    Stats.MLT_STARTED_DISPLAY_UPDATE();

    camera.film.updateDisplay(x0, y0, x1, y1, splatScale);
    camera.film.writeImage(splatScale);

    Stats.MLT_FINISHED_DISPLAY_UPDATE();
    Stats.MLT_FINISHED_MLT_TASK(this);
  }

  int taskNum;
  double dx, dy;
  int currentPixelSample;
  int x0, x1, y0, y1;
  double t0, t1;
  double b;
  _MLTSample initialSample;
  Scene scene;
  Camera camera;
  MetropolisRenderer renderer;
  Distribution1D lightDistribution;
}
