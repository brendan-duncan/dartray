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
part of dartray;

/**
 * Manages the rendering process, either rendering locally or submitting the
 * job to one or more isolates (web workers).
 */
abstract class RenderManagerInterface extends ResourceManager {
  DartRay dartray;
  RenderIsolate isolate;
  String scenePath;
  OutputImage renderOutput;

  RenderManagerInterface(this.scenePath) {
    RegisterStandardPlugins();
    dartray = new DartRay(this);
  }

  static void RegisterStandardPlugins() {
    // If 'sphere' has been registered, we can assume the rest of the standard
    // plugins have been registered too.
    if (Plugin.shape('sphere') != null) {
      return;
    }

    Plugin.registerAccelerator('bvh', BVHAccel.Create);
    Plugin.registerAccelerator('grid', GridAccel.Create);
    Plugin.registerAccelerator('kdtree', KdTreeAccel.Create);
    Plugin.registerAccelerator('bruteforce', BruteForceAccel.Create);

    Plugin.registerCamera('environment', EnvironmentCamera.Create);
    Plugin.registerCamera('orthographic', OrthographicCamera.Create);
    Plugin.registerCamera('perspective', PerspectiveCamera.Create);

    Plugin.registerFilm('image', ImageFilm.Create);

    Plugin.registerFilter('box', BoxFilter.Create);
    Plugin.registerFilter('gaussian', GaussianFilter.Create);
    Plugin.registerFilter('sinc', LanczosSincFilter.Create);
    Plugin.registerFilter('mitchell', MitchellFilter.Create);
    Plugin.registerFilter('triangle', TriangleFilter.Create);

    Plugin.registerSurfaceIntegrator('ambientocclusion',
        AmbientOcclusionIntegrator.Create);
    Plugin.registerSurfaceIntegrator('diffuseprt', DiffusePRTIntegrator.Create);
    Plugin.registerSurfaceIntegrator('directlighting',
        DirectLightingIntegrator.Create);
    Plugin.registerSurfaceIntegrator('glossyprt', GlossyPRTIntegrator.Create);
    Plugin.registerSurfaceIntegrator('igi', IGIIntegrator.Create);
    Plugin.registerSurfaceIntegrator('irradiancecache',
        IrradianceCacheIntegrator.Create);
    Plugin.registerSurfaceIntegrator('path', PathIntegrator.Create);
    Plugin.registerSurfaceIntegrator('photonmap', PhotonMapIntegrator.Create);
    Plugin.registerSurfaceIntegrator('exphotonmap', PhotonMapIntegrator.Create);
    Plugin.registerSurfaceIntegrator('whitted', WhittedIntegrator.Create);
    Plugin.registerSurfaceIntegrator('useprobes', UseProbesIntegrator.Create);
    Plugin.registerSurfaceIntegrator('dipolesubsurface',
        DipoleSubsurfaceIntegrator.Create);

    Plugin.registerLight('distant', DistantLight.Create);
    Plugin.registerLight('point', PointLight.Create);
    Plugin.registerLight('spot', SpotLight.Create);
    Plugin.registerLight('infinite', InfiniteAreaLight.Create);
    Plugin.registerLight('goniometric', GoniometricLight.Create);
    Plugin.registerLight('projection', ProjectionLight.Create);

    Plugin.registerAreaLight('diffuse', DiffuseAreaLight.Create);
    Plugin.registerAreaLight('area', DiffuseAreaLight.Create);

    Plugin.registerMaterial('glass', GlassMaterial.Create);
    Plugin.registerMaterial('kdsubsurface', KdSubsurfaceMaterial.Create);
    Plugin.registerMaterial('matte', MatteMaterial.Create);
    Plugin.registerMaterial('measured', MeasuredMaterial.Create);
    Plugin.registerMaterial('metal', MetalMaterial.Create);
    Plugin.registerMaterial('mirror', MirrorMaterial.Create);
    Plugin.registerMaterial('plastic', PlasticMaterial.Create);
    Plugin.registerMaterial('shinymetal', ShinyMetalMaterial.Create);
    Plugin.registerMaterial('substrate', SubstrateMaterial.Create);
    Plugin.registerMaterial('subsurface', SubsurfaceMaterial.Create);
    Plugin.registerMaterial('translucent', TranslucentMaterial.Create);
    Plugin.registerMaterial('uber', UberMaterial.Create);

    Plugin.registerPixelSampler('linear', LinearPixelSampler.Create);
    Plugin.registerPixelSampler('random', RandomPixelSampler.Create);
    Plugin.registerPixelSampler('tile', TilePixelSampler.Create);

    Plugin.registerSampler('adaptive', AdaptiveSampler.Create);
    Plugin.registerSampler('bestcandidate', BestCandidateSampler.Create);
    Plugin.registerSampler('halton', HaltonSampler.Create);
    Plugin.registerSampler('lowdiscrepancy', LowDiscrepancySampler.Create);
    Plugin.registerSampler('random', RandomSampler.Create);
    Plugin.registerSampler('stratified', StratifiedSampler.Create);

    Plugin.registerShape('cone', Cone.Create);
    Plugin.registerShape('cylinder', Cylinder.Create);
    Plugin.registerShape('disk', Disk.Create);
    Plugin.registerShape('heightfield', Heightfield.Create);
    Plugin.registerShape('hyperboloid', Hyperboloid.Create);
    Plugin.registerShape('loopsubdiv', LoopSubdivision.Create);
    Plugin.registerShape('nurbs', Nurbs.Create);
    Plugin.registerShape('paraboloid', Paraboloid.Create);
    Plugin.registerShape('sphere', Sphere.Create);
    Plugin.registerShape('trianglemesh', TriangleMesh.Create);

    Plugin.registerFloatTexture('bilerp', BilerpTexture.CreateFloat);
    Plugin.registerSpectrumTexture('bilerp', BilerpTexture.CreateSpectrum);
    Plugin.registerFloatTexture('checkerboard',
        CheckerboardTexture.CreateFloat);
    Plugin.registerSpectrumTexture('checkerboard',
        CheckerboardTexture.CreateSpectrum);
    Plugin.registerFloatTexture('constant', ConstantTexture.CreateFloat);
    Plugin.registerSpectrumTexture('constant', ConstantTexture.CreateSpectrum);
    Plugin.registerFloatTexture('dots', DotsTexture.CreateFloat);
    Plugin.registerSpectrumTexture('dots', DotsTexture.CreateSpectrum);
    Plugin.registerFloatTexture('fbm', FBmTexture.CreateFloat);
    Plugin.registerSpectrumTexture('fbm', FBmTexture.CreateSpectrum);
    Plugin.registerFloatTexture('imagemap', ImageTexture.CreateFloat);
    Plugin.registerSpectrumTexture('imagemap', ImageTexture.CreateSpectrum);
    Plugin.registerFloatTexture('marble', MarbleTexture.CreateFloat);
    Plugin.registerSpectrumTexture('marble', MarbleTexture.CreateSpectrum);
    Plugin.registerFloatTexture('mix', MixTexture.CreateFloat);
    Plugin.registerSpectrumTexture('mix', MixTexture.CreateSpectrum);
    Plugin.registerFloatTexture('scale', ScaleTexture.CreateFloat);
    Plugin.registerSpectrumTexture('scale', ScaleTexture.CreateSpectrum);
    Plugin.registerFloatTexture('uv', UVTexture.CreateFloat);
    Plugin.registerSpectrumTexture('uv', UVTexture.CreateSpectrum);
    Plugin.registerFloatTexture('windy', WindyTexture.CreateFloat);
    Plugin.registerSpectrumTexture('windy', WindyTexture.CreateSpectrum);
    Plugin.registerFloatTexture('wrinkled', WrinkledTexture.CreateFloat);
    Plugin.registerSpectrumTexture('wrinkled', WrinkledTexture.CreateSpectrum);

    Plugin.registerVolumeIntegrator('emission', EmissionIntegrator.Create);
    Plugin.registerVolumeIntegrator('single',
        SingleScatteringIntegrator.Create);

    Plugin.registerVolumeRegion('exponential', ExponentialDensityRegion.Create);
    Plugin.registerVolumeRegion('homogeneous', HomogeneousVolumeRegion.Create);
    Plugin.registerVolumeRegion('volumegrid', VolumeGridDensity.Create);
  }

  /**
   * This is called from an Isolate to initialize the RenderManager for the
   * Isolate.
   */
  void startIsolate([SendPort port]) {
    if (port != null) {
      isolate = new RenderIsolate(this);
      isolate.start(port);
    }
  }

  void pause() {
    if (isolates == null) {
      return;
    }

    for (RenderTask task in isolates) {
      task.pause();
    }
  }

  void resume() {
    if (isolates == null) {
      return;
    }

    for (RenderTask task in isolates) {
      task.resume;
    }
  }

  /**
   * [isolateUri] is the path to the script that initiates the
   * RenderIsolate.  This script usually would do nothing more than:
   * import 'dart:isolate';
   * import 'package:dartray/dartray.dart';
   * void main(List<String> args, SendPort port) {
   *   new RenderIsolate().start(port);
   * }
   */
  Future<OutputImage> render(String path, {String isolate,
              LogCallback log, PreviewCallback preview,
              RenderOverrides overrides, int numThreads: 1}) {
    if (log != null) {
      Log = log;
    }

    if (path.contains('/')) {
      int i = path.lastIndexOf('/');
      scenePath = path.substring(0, i);
      path = path.substring(i + 1);
    }

    Completer<OutputImage> completer = new Completer<OutputImage>();

    if (isolate == null) {
      LogInfo('STARTING DartRay Render');
      dartray.renderScene(path, overrides: overrides).then((output) {
        completer.complete(output);
      });

      return completer.future;
    }

    int tasksRemaining = numThreads;
    isolates = new List<RenderTask>(numThreads);

    for (int i = 0; i < numThreads; ++i) {
      isolates[i] = new RenderTask(preview, i, numThreads);
      isolates[i].render(path, isolate, overrides: overrides).then((output) {
        if (numThreads > 1) {
          if (renderOutput == null ||
              renderOutput.imageWidth != output.imageWidth ||
              renderOutput.imageHeight != output.imageHeight) {
            renderOutput = new OutputImage(0, 0, output.imageWidth,
                                           output.imageHeight);
          }

          for (int y = 0; y < output.height; ++y) {
            int pi = ((output.yOffset + y) * renderOutput.imageWidth * 3) +
                     (output.xOffset * 3);
            for (int x = 0; x < output.width; ++x, pi += 3) {
              renderOutput.rgb[pi] = output.rgb[pi];
              renderOutput.rgb[pi + 1] = output.rgb[pi + 1];
              renderOutput.rgb[pi + 2] = output.rgb[pi + 2];
            }
          }
        } else {
          renderOutput = output;
        }

        tasksRemaining--;
        if (tasksRemaining == 0) {
          completer.complete(renderOutput);
        }
      }, onError: (msg) {
        LogError('ERROR Thread $i: $msg');
        tasksRemaining--;
        /*if (!completer.isCompleted) {
          completer.complete(output);
        }*/
      });
    }

    return completer.future;
  }

  List<RenderTask> isolates;
}

