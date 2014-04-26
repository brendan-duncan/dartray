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
library dartray;

import 'dart:async';
import 'dart:isolate';
import 'dart:math' as Math;
import 'dart:typed_data';

import 'package:archive/archive.dart';
import 'package:image/image.dart';

import '../accelerators/accelerators.dart';
import '../core/core.dart';
import '../cameras/cameras.dart';
import '../film/film.dart';
import '../filters/filters.dart';
import '../lights/lights.dart';
import '../materials/materials.dart';
import '../pixel_samplers/pixel_samplers.dart';
import '../renderers/renderers.dart';
import '../samplers/samplers.dart';
import '../shapes/shapes.dart';
import '../surface_integrators/surface_integrators.dart';
import '../textures/textures.dart';
import '../volume_integrators/volume_integrators.dart';
import '../volume_regions/volume_regions.dart';

part 'graphics_state.dart';
part 'pbrt_lexer.dart';
part 'pbrt_parser.dart';
part 'render_isolate.dart';
part 'render_manager_interface.dart';
part 'render_options.dart';
part 'render_task.dart';
part 'transform_set.dart';


/**
 * [DartRay] provides an API and file format compatible with the PBRT
 * rendering system API and file format.
 */
class DartRay {
  ResourceManager resourceManager;
  RenderOverrides overrides;
  OutputImage outputImage;

  DartRay(this.resourceManager);

  Future<OutputImage> renderScene(String scene, {RenderOverrides overrides}) {
    Stopwatch t = new Stopwatch()..start();
    Completer<OutputImage> c = new Completer<OutputImage>();

    this.overrides = overrides;

    if (this.overrides != null) {
      LogDebug('OVERRIDES: ${this.overrides.toJson()}');
    }

    try {
      PbrtParser parser = new PbrtParser(this);
      parser.parse(scene).then((x) {
        c.complete(outputImage);
      });
    } catch (e) {
      LogError('EXCEPTION: $e');
      c.completeError(e);
    }

    return c.future;
  }

  static const int MAX_TRANSFORMS = 2;
  static const int START_TRANSFORM_BITS = (1 << 0);
  static const int END_TRANSFORM_BITS = (1 << 1);
  static const int ALL_TRANSFORMS_BITS = ((1 << MAX_TRANSFORMS) - 1);

  static const int STATE_UNINITIALIZED = 0;
  static const int STATE_OPTIONS_BLOCK = 1;
  static const int STATE_WORLD_BLOCK = 2;

  int _currentApiState = STATE_OPTIONS_BLOCK;
  TransformSet _curTransform = new TransformSet();
  int _activeTransformBits = ALL_TRANSFORMS_BITS;
  Map<String, TransformSet> _namedCoordinateSystems = {};
  RenderOptions _renderOptions = new RenderOptions();
  GraphicsState _graphicsState = new GraphicsState();
  List<GraphicsState> _pushedGraphicsStates = [];
  List<TransformSet> _pushedTransforms = [];
  List<int> _pushedActiveTransformBits = [];

  // API Cleanup
  void cleanup() {
    //ProbesCleanup();
    if (_currentApiState == STATE_UNINITIALIZED) {
      LogWarning('cleanup() called without init().');
    } else if (_currentApiState == STATE_WORLD_BLOCK) {
      LogWarning('cleanup() called while inside world block.');
    }

    _currentApiState = STATE_UNINITIALIZED;
  }

  void setTask(int taskNum, int taskCount) {
    _renderOptions.taskNum = taskNum;
    _renderOptions.taskCount = taskCount;
  }

  void identity() {
    for (int i = 0; i < MAX_TRANSFORMS; ++i) {
      if (_activeTransformBits & (1 << i) != 0) {
        _curTransform[i] = new Transform();
      }
    }
  }

  void translate(double dx, double dy, double dz) {
    for (int i = 0; i < MAX_TRANSFORMS; ++i) {
      if (_activeTransformBits & (1 << i) != 0) {
        _curTransform[i] = _curTransform[i] *
                          Transform.Translate(new Vector(dx, dy, dz));
      }
    }
  }

  void transform(Matrix4x4 tr) {
    for (int i = 0; i < MAX_TRANSFORMS; ++i) {
      if (_activeTransformBits & (1 << i) != 0) {
        _curTransform[i] = new Transform(tr);
      }
    }
  }

  void concatTransform(Matrix4x4 tr) {
    for (int i = 0; i < MAX_TRANSFORMS; ++i) {
      if (_activeTransformBits & (1 << i) != 0) {
        _curTransform[i] = _curTransform[i] * new Transform(tr);
      }
    }
  }

  void rotate(double angle, double dx, double dy, double dz) {
    for (int i = 0; i < MAX_TRANSFORMS; ++i) {
      if (_activeTransformBits & (1 << i) != 0) {
        _curTransform[i] = _curTransform[i] *
                          Transform.Rotate(angle, new Vector(dx, dy, dz));
      }
    }
  }

  void scale(double sx, double sy, double sz) {
    for (int i = 0; i < MAX_TRANSFORMS; ++i) {
      if (_activeTransformBits & (1 << i) != 0) {
        _curTransform[i] = _curTransform[i] * Transform.Scale(sx, sy, sz);
      }
    }
  }

  void lookAt(double ex, double ey, double ez, double lx, double ly,
              double lz, double ux, double uy, double uz) {
    for (int i = 0; i < MAX_TRANSFORMS; ++i) {
      if (_activeTransformBits & (1 << i) != 0) {
        _curTransform[i] = _curTransform[i] *
            Transform.LookAt(new Point(ex, ey, ez), new Point(lx, ly, lz),
                             new Vector(ux, uy, uz));
      }
    }
  }

  void coordinateSystem(String name) {
    _namedCoordinateSystems[name] = _curTransform;
  }

  void coordSysTransform(String name) {
    if (_namedCoordinateSystems.containsKey(name)) {
      _curTransform = new TransformSet.from(_namedCoordinateSystems[name]);
    } else {
      LogWarning('Couldn\'t find named coordinate system \'$name\'');
    }
  }

  void activeTransformAll() {
    _activeTransformBits = ALL_TRANSFORMS_BITS;
  }

  void activeTransformEndTime() {
    _activeTransformBits = END_TRANSFORM_BITS;
  }

  void activeTransformStartTime() {
    _activeTransformBits = START_TRANSFORM_BITS;
  }

  void transformTimes(double start, double end) {
    _renderOptions.transformStartTime = start;
    _renderOptions.transformEndTime = end;
  }

  void pixelFilter(String name, ParamSet params) {
    _renderOptions.filterName = name;
    _renderOptions.filterParams = params;
  }

  void film(String type, ParamSet params) {
    _renderOptions.filmParams = params;
    _renderOptions.filmName = type;
  }

  void setPreviewCallback(PreviewCallback cb) {
    _renderOptions.previewCallback = cb;
  }

  void pixels(String name, ParamSet params) {
    _renderOptions.pixelSamplerName = name;
    _renderOptions.pixelSamplerParams = params;
  }

  void sampler(String name, ParamSet params) {
    _renderOptions.samplerName = name;
    _renderOptions.samplerParams = params;
  }

  void accelerator(String name, ParamSet params) {
    _renderOptions.acceleratorName = name;
    _renderOptions.acceleratorParams = params;
  }

  void surfaceIntegrator(String name, ParamSet params) {
    _renderOptions.surfaceIntegratorName = name;
    _renderOptions.surfaceIntegratorParams = params;
  }

  void volumeIntegrator(String name, ParamSet params) {
    _renderOptions.volumeIntegratorName = name;
    _renderOptions.volumeIntegratorParams = params;
  }

  void renderer(String name, ParamSet params) {
    _renderOptions.rendererName = name;
    _renderOptions.rendererParams = params;
  }

  void camera(String name, ParamSet params) {
    _renderOptions.cameraName = name;
    _renderOptions.cameraParams = params;
    _renderOptions.cameraToWorld = TransformSet.Inverse(_curTransform);
    _namedCoordinateSystems['camera'] = _renderOptions.cameraToWorld;
  }

  void worldBegin() {
    _currentApiState = STATE_WORLD_BLOCK;
    for (int i = 0; i < MAX_TRANSFORMS; ++i) {
      _curTransform[i] = new Transform();
    }
    _activeTransformBits = ALL_TRANSFORMS_BITS;
    _namedCoordinateSystems['world'] = _curTransform;
  }

  void attributeBegin() {
    _pushedGraphicsStates.add(new GraphicsState.from(_graphicsState));
    _pushedTransforms.add(new TransformSet.from(_curTransform));
    _pushedActiveTransformBits.add(_activeTransformBits);
  }

  void attributeEnd() {
    if (_pushedGraphicsStates.isEmpty) {
      LogWarning('Unmatched attributeEnd() encountered. '
                 'Ignoring it.');
      return;
    }
    _graphicsState = _pushedGraphicsStates.removeLast();
    _curTransform = _pushedTransforms.removeLast();
    _activeTransformBits = _pushedActiveTransformBits.removeLast();
  }

  void transformBegin() {
    _pushedTransforms.add(new TransformSet.from(_curTransform));
    _pushedActiveTransformBits.add(_activeTransformBits);
  }

  void transformEnd() {
    if (_pushedTransforms.isEmpty) {
      LogWarning('Unmatched pbrtTransformEnd() encountered. '
                 'Ignoring it.');
      return;
    }

    _curTransform = _pushedTransforms.removeLast();
    _activeTransformBits = _pushedActiveTransformBits.removeLast();
  }

  void texture(String name, String type, String texname, ParamSet params) {
    TextureParams tp = new TextureParams(params, params,
                            _graphicsState.doubleTextures,
                            _graphicsState.spectrumTextures);
    if (type == 'float')  {
      // Create _double_ texture and store in _doubleTextures_
      if (_graphicsState.doubleTextures.containsKey(name)) {
        LogWarning('Texture \'$name\' being redefined');
      }

      //WARN_IF_ANIMATED_TRANSFORM('Texture');
      Texture ft = _makeFloatTexture(texname, _curTransform[0], tp);
      if (ft != null) {
        _graphicsState.doubleTextures[name] = ft;
      }
    } else if (type == 'color' || type == 'spectrum')  {
      // Create _color_ texture and store in _spectrumTextures_
      if (_graphicsState.spectrumTextures.containsKey(name)) {
        LogWarning('Texture \'$name\' being redefined');
      }

      //WARN_IF_ANIMATED_TRANSFORM('Texture');
      Texture st = _makeSpectrumTexture(texname, _curTransform[0], tp);
      if (st != null) {
        _graphicsState.spectrumTextures[name] = st;
      }
    } else {
      LogWarning('Texture type \'$type\' unknown.');
    }
  }

  void material(String name, ParamSet params) {
    _graphicsState.material = name;
    _graphicsState.materialParams = params;
    _graphicsState.currentNamedMaterial = '';
  }

  void makeNamedMaterial(String name, ParamSet params) {
    // error checking, warning if replace, what to use for transform?
    TextureParams mp = new TextureParams(params, _graphicsState.materialParams,
                     _graphicsState.doubleTextures,
                     _graphicsState.spectrumTextures);
    String matName = mp.findString('type');

    //WARN_IF_ANIMATED_TRANSFORM('MakeNamedMaterial');
    if (matName.isEmpty) {
      LogWarning('No parameter string \'type\' found in MakeNamedMaterial');
    } else {
      Material mtl = _makeMaterial(matName, _curTransform[0], mp);
      if (mtl != null) {
        _graphicsState.namedMaterials[name] = mtl;
      }
    }
  }

  void namedMaterial(String name) {
    _graphicsState.currentNamedMaterial = name;
  }

  void lightSource(String name, ParamSet params) {
    //WARN_IF_ANIMATED_TRANSFORM('LightSource');
    Light lt = _makeLight(name, _curTransform[0], params);
    if (lt == null) {
      LogWarning('pbrtLightSource: light type \'$name\' unknown.');
    } else {
      _renderOptions.lights.add(lt);
    }
  }

  void areaLightSource(String name, ParamSet params) {
    _graphicsState.areaLight = name;
    _graphicsState.areaLightParams = params;
  }

  void shape(String name, ParamSet params) {
    Primitive prim;
    AreaLight area;

    if (!_curTransform.isAnimated()) {
      // Create primitive for static shape
      Transform obj2world = new Transform.from(_curTransform[0]);
      Transform world2obj = Transform.Inverse(obj2world);

      Shape shape = _makeShape(name, obj2world, world2obj,
                              _graphicsState.reverseOrientation, params);
      if (shape == null) {
        return;
      }

      Material mtl = _createMaterial(params);

      params.reportUnused();

      // Possibly create area light for shape
      if (_graphicsState.areaLight != '') {
        area = _makeAreaLight(_graphicsState.areaLight, _curTransform[0],
                             _graphicsState.areaLightParams, shape);
      }

      prim = new GeometricPrimitive(shape, mtl, area);
    } else {
      // Create primitive for animated shape

      // Create initial _Shape_ for animated shape
      if (_graphicsState.areaLight != '') {
        LogWarning('Ignoring currently set area light when creating '
              'animated shape');
      }

      Transform identity = new Transform();

      Shape shape = _makeShape(name, identity, identity,
                              _graphicsState.reverseOrientation, params);
      if (shape == null) {
        return;
      }

      Material mtl = _createMaterial(params);

      params.reportUnused();

      // Get _animatedWorldToObject_ transform for shape
      assert(MAX_TRANSFORMS == 2);

      Transform world2obj0 = Transform.Inverse(_curTransform[0]);
      Transform world2obj1 = Transform.Inverse(_curTransform[1]);

      AnimatedTransform animatedWorldToObject =
          new AnimatedTransform(world2obj0, _renderOptions.transformStartTime,
                                world2obj1, _renderOptions.transformEndTime);

      Primitive baseprim = new GeometricPrimitive(shape, mtl, null);

      if (!baseprim.canIntersect()) {
        // Refine animated shape and create BVH if more than one shape created
        List<Primitive> refinedPrimitives = [];
        baseprim.fullyRefine(refinedPrimitives);
        if (refinedPrimitives.isEmpty) {
          return;
        }
        if (refinedPrimitives.length > 1) {
          baseprim = new BVHAccel(refinedPrimitives);
        } else {
          baseprim = refinedPrimitives[0];
        }
      }

      prim = new TransformedPrimitive(baseprim, animatedWorldToObject);
    }

    // Add primitive to scene or current instance
    if (_renderOptions.currentInstance != null) {
      if (area != null) {
        LogWarning('Area lights not supported with object instancing');
      }

      _renderOptions.currentInstance.add(prim);
    } else {
      _renderOptions.primitives.add(prim);
      if (area != null) {
        _renderOptions.lights.add(area);
      }
    }
  }

  void reverseOrientation() {
    _graphicsState.reverseOrientation = !_graphicsState.reverseOrientation;
  }

  void volume(String name, ParamSet params) {
    VolumeRegion vr = _makeVolumeRegion(name, _curTransform[0], params);
    if (vr != null) {
      _renderOptions.volumeRegions.add(vr);
    }
  }

  void objectBegin(String name) {
    attributeBegin();
    if (_renderOptions.currentInstance != null) {
      LogWarning('ObjectBegin called inside of instance definition');
    }
    _renderOptions.instances[name] = new List<Primitive>();
    _renderOptions.currentInstance = _renderOptions.instances[name];
  }


  void objectEnd() {
    if (_renderOptions.currentInstance == null) {
      LogWarning('ObjectEnd called outside of instance definition');
    }
    _renderOptions.currentInstance = null;
    attributeEnd();
  }


  void objectInstance(String name) {
    // Object instance error checking
    if (_renderOptions.currentInstance != null) {
      LogWarning('ObjectInstance can\'t be called inside instance definition');
      return;
    }

    if (!_renderOptions.instances.containsKey(name)) {
      LogWarning('Unable to find instance named \'$name\'');
      return;
    }

    List<Primitive> inst = _renderOptions.instances[name];

    if (inst.isEmpty) {
      return;
    }

    if (inst.length > 1 || !inst[0].canIntersect()) {
      // Refine instance _Primitive_s and create aggregate
      Primitive accel = _makeAccelerator(_renderOptions.acceleratorName,
                                         inst,
                                         _renderOptions.acceleratorParams);

      if (accel == null) {
        accel = _makeAccelerator('bvh', inst, new ParamSet());
      }

      if (accel == null) {
        LogSevere('Unable to create \'bvh\' accelerator');
      }

      inst.clear();
      inst.add(accel);
    }

    assert(MAX_TRANSFORMS == 2);

    Transform world2instance0 = Transform.Inverse(_curTransform[0]);
    Transform world2instance1 = Transform.Inverse(_curTransform[1]);

    AnimatedTransform animatedWorldToInstance = new AnimatedTransform(
        world2instance0, _renderOptions.transformStartTime,
        world2instance1, _renderOptions.transformEndTime);

    Primitive prim = new TransformedPrimitive(inst[0], animatedWorldToInstance);
    _renderOptions.primitives.add(prim);
  }


  Future worldEnd() {
    // Ensure there are no pushed graphics states
    while (_pushedGraphicsStates.isNotEmpty) {
      LogWarning('Missing end to pbrtAttributeBegin()');
      _pushedGraphicsStates.removeLast();
      _pushedTransforms.removeLast();
    }

    while (_pushedTransforms.isNotEmpty) {
      LogWarning('Missing end to pbrtTransformBegin()');
      _pushedTransforms.removeLast();
    }

    // Create scene and render
    Renderer renderer = _makeRenderer();
    Scene scene = _makeScene();

    Future future;

    if (renderer != null && scene != null) {
      Completer c = new Completer();
      future = c.future;

      resourceManager.waitUntilReady().then((_) {
        try {
          future = renderer.render(scene);
          future.then((OutputImage output) {
            outputImage = output;
            c.complete();
          }).catchError((e) {
            c.completeError(e);
          });
        } catch (e) {
          c.completeError(e);
        }
      });
    }

    // Clean up after rendering
    _graphicsState = new GraphicsState();

    _currentApiState = STATE_OPTIONS_BLOCK;
    //ProbesPrint(stdout);
    for (int i = 0; i < MAX_TRANSFORMS; ++i) {
      _curTransform.t[i] = new Transform();
    }

    _activeTransformBits = ALL_TRANSFORMS_BITS;
    _namedCoordinateSystems.clear();

    return future;
  }


  Scene _makeScene() {
    // Initialize _volumeRegion_ from volume region(s)
    VolumeRegion volumeRegion;
    if (_renderOptions.volumeRegions.length == 0) {
      volumeRegion = null;
    } else if (_renderOptions.volumeRegions.length == 1) {
      volumeRegion = _renderOptions.volumeRegions[0];
    } else {
      volumeRegion = new AggregateVolume(_renderOptions.volumeRegions);
    }

    Primitive accelerator = _makeAccelerator(_renderOptions.acceleratorName,
                                             _renderOptions.primitives,
                                             _renderOptions.acceleratorParams);

    if (accelerator == null) {
      accelerator = _makeAccelerator('bvh', _renderOptions.primitives,
                                     new ParamSet());
    }

    if (accelerator == null) {
      LogSevere('Unable to create \'bvh\' accelerator.');
    }

    Scene scene = new Scene(accelerator, _renderOptions.lights, volumeRegion);

    // Erase primitives, lights, and volume regions from _RenderOptions_
    _renderOptions.primitives.clear();
    _renderOptions.lights.clear();
    _renderOptions.volumeRegions.clear();

    return scene;
  }

  Renderer _makeRenderer() {
    Renderer renderer;

    Filter filter = _makeFilter(_renderOptions.filterName,
                                _renderOptions.filterParams);

    Film film = _makeFilm(_renderOptions.filmName,
                          _renderOptions.filmParams,
                          filter,
                          _renderOptions.previewCallback);
    if (film == null) {
      LogSevere('Unable to create film.');
    }

    Camera camera = _makeCamera(_renderOptions.cameraName,
                                _renderOptions.cameraParams,
                                _renderOptions.cameraToWorld,
                                _renderOptions.transformStartTime,
                                _renderOptions.transformEndTime,
                                film);

    if (camera == null) {
      LogSevere('Unable to create camera.');
    }

    String name = _renderOptions.rendererName;
    ParamSet paramSet = _renderOptions.rendererParams;

    if (overrides != null && overrides.rendererName != null) {
      name = overrides.rendererName;
      paramSet = overrides.rendererParams;
    }

    if (name == 'createprobes') {
      // Create surface and volume integrators
      SurfaceIntegrator surfaceIntegrator =
          _makeSurfaceIntegrator(_renderOptions.surfaceIntegratorName,
                                 _renderOptions.surfaceIntegratorParams);
      if (surfaceIntegrator == null) {
        LogSevere('Unable to create surface integrator.');
      }

      VolumeIntegrator volumeIntegrator =
          _makeVolumeIntegrator(_renderOptions.volumeIntegratorName,
                                _renderOptions.volumeIntegratorParams);
      if (volumeIntegrator == null) {
        LogSevere('Unable to create volume integrator.');
      }

      renderer = CreateProbesRenderer.Create(camera, surfaceIntegrator,
                                             volumeIntegrator, paramSet);

      paramSet.reportUnused();

      // Warn if no light sources are defined
      if (_renderOptions.lights.isEmpty) {
        LogWarning('No light sources defined in scene; '
                   'possibly rendering a black image.');
      }
    } else if (name == 'surfacepoints') {
      Point pCamera = camera.cameraToWorld.transformPoint(camera.shutterOpen,
                                                          Point.ZERO);

      renderer = SurfacePointsRenderer.Create(paramSet, pCamera,
                                              camera.shutterOpen);

      paramSet.reportUnused();
    } else if (name == 'aggregatetest') {
      renderer = AggregateTestRenderer.Create(paramSet,
                                              _renderOptions.primitives);
      paramSet.reportUnused();
    } else if (name == 'metropolis') {
      renderer = MetropolisRenderer.Create(paramSet,
                                           camera, _renderOptions.taskNum,
                                           _renderOptions.taskCount);
      paramSet.reportUnused();

      // Warn if no light sources are defined
      if (_renderOptions.lights.length == 0) {
        LogWarning('No light sources defined in scene; '
              'possibly rendering a black image.');
      }
    } else {
      if (name != 'sampler') {
        LogWarning('Renderer type \'${name}\' unknown. '
                   'Using \'sampler\'.');
      }

      paramSet.reportUnused();

      PixelSampler pixels = _makePixelSampler(_renderOptions.pixelSamplerName,
                                           _renderOptions.pixelSamplerParams);

      Sampler sampler = _makeSampler(_renderOptions.samplerName,
                                     _renderOptions.samplerParams,
                                     camera.film, camera, pixels);
      if (sampler == null) {
        LogSevere('Unable to create sampler.');
      }

      // Create surface and volume integrators
      SurfaceIntegrator surfaceIntegrator =
          _makeSurfaceIntegrator(_renderOptions.surfaceIntegratorName,
                                 _renderOptions.surfaceIntegratorParams);
      if (surfaceIntegrator == null) {
        LogSevere('Unable to create surface integrator.');
      }

      VolumeIntegrator volumeIntegrator =
          _makeVolumeIntegrator(_renderOptions.volumeIntegratorName,
                                _renderOptions.volumeIntegratorParams);
      if (volumeIntegrator == null) {
        LogSevere('Unable to create volume integrator.');
      }

      renderer = new SamplerRenderer(sampler, camera, surfaceIntegrator,
                                     volumeIntegrator, _renderOptions.taskNum,
                                     _renderOptions.taskCount);

      // Warn if no light sources are defined
      if (_renderOptions.lights.isEmpty) {
        LogWarning('No light sources defined in scene; '
              'possibly rendering a black image.');
      }
    }

    return renderer;
  }

  Shape _makeShape(String name, Transform object2world, Transform world2object,
                  bool reverseOrientation, ParamSet paramSet) {
    ShapeCreator sh = Plugin.shape(name);
    if (sh == null) {
      LogWarning('Shape \'${name}\' unknown.');
      return null;
    }

    Shape s = sh(object2world, world2object, reverseOrientation, paramSet);
    paramSet.reportUnused();

    return s;
  }

  Material _createMaterial(ParamSet params) {
    TextureParams mp = new TextureParams(params, _graphicsState.materialParams,
        _graphicsState.doubleTextures,
        _graphicsState.spectrumTextures);

    Material mtl;
    if (_graphicsState.currentNamedMaterial != '' &&
        _graphicsState.namedMaterials.containsKey(_graphicsState.currentNamedMaterial)) {
      mtl = _graphicsState.namedMaterials[_graphicsState.currentNamedMaterial];
    }

    if (mtl == null) {
      mtl = _makeMaterial(_graphicsState.material, _curTransform[0], mp);
    }

    if (mtl == null) {
      mtl = _makeMaterial('matte', _curTransform[0], mp);
    }

    if (mtl == null) {
      LogSevere('Unable to create \'matte\' material?!');
    }

    return mtl;
  }

  Material _makeMaterial(String name, Transform mtl2world, TextureParams mp) {
    // Special case the 'mix' material
    if (name == 'mix') {
      String m1 = mp.findString('namedmaterial1', '');
      String m2 = mp.findString('namedmaterial2', '');
      Material mat1 = _graphicsState.namedMaterials[m1];
      Material mat2 = _graphicsState.namedMaterials[m2];
      if (mat1 == null) {
        LogWarning('Named material \'$m1\' undefined.  Using \'matte\'');
        mat1 = _makeMaterial('matte', _curTransform[0], mp);
      }
      if (mat2 == null) {
        LogWarning('Named material \'$m2\' undefined.  Using \'matte\'');
        mat2 = _makeMaterial('matte', _curTransform[0], mp);
      }

      Material m = MixMaterial.Create(mtl2world, mp, mat1, mat2);
      mp.reportUnused();

      return m;
    }

    if (Plugin.material(name) == null) {
      LogWarning('Material \'$name\' unknown.');
      return null;
    }

    Material material = Plugin.material(name)(mtl2world, mp);
    mp.reportUnused();

    return material;
  }

  Texture _makeFloatTexture(String name, Transform tex2world,
                            TextureParams tp) {
    if (Plugin.floatTexture(name) == null) {
      LogWarning('Texture \'${name}\' unknown.');
      return ConstantTexture.CreateFloat(tex2world, tp);
    }

    Texture t = Plugin.floatTexture(name)(tex2world, tp);
    tp.reportUnused();

    return t;
  }

  Texture _makeSpectrumTexture(String name, Transform tex2world,
                              TextureParams tp) {
    if (Plugin.spectrumTexture(name) == null) {
      LogWarning('Texture \'${name}\' unknown.');
      return ConstantTexture.CreateSpectrum(tex2world, tp);
    }

    Texture t = Plugin.spectrumTexture(name)(tex2world, tp);
    tp.reportUnused();

    return t;
  }

  Light _makeLight(String name, Transform light2world, ParamSet paramSet) {
    if (Plugin.light(name) == null) {
      LogWarning('Light \'${name}\' unknown.');
      return null;
    }

    Light l = Plugin.light(name)(light2world, paramSet);
    paramSet.reportUnused();

    return l;
  }

  AreaLight _makeAreaLight(String name, Transform light2world,
                           ParamSet paramSet, Shape shape) {
    if (Plugin.areaLight(name) == null) {
      LogWarning('Area Light \'${name}\' unknown.');
      return null;
    }

    Light l = Plugin.areaLight(name)(light2world, paramSet, shape);
    paramSet.reportUnused();

    return l;
  }

  VolumeRegion _makeVolumeRegion(String name, Transform volume2world,
                                ParamSet paramSet) {
    if (Plugin.volumeRegion(name) == null) {
      LogWarning('Volume Region \'${name}\' unknown.');
      return null;
    }

    VolumeRegion v = Plugin.volumeRegion(name)(volume2world, paramSet);
    paramSet.reportUnused();

    return v;
  }

  SurfaceIntegrator _makeSurfaceIntegrator(String name, ParamSet paramSet) {
    if (overrides != null && overrides.volumeIntegratorName != null) {
      name = overrides.surfaceIntegratorName;
      paramSet = overrides.surfaceIntegratorParams;
    }

    if (Plugin.surfaceIntegrator(name) == null) {
      LogWarning('Surface Integrator \'${name}\' unknown.');
      return null;
    }

    SurfaceIntegrator si = Plugin.surfaceIntegrator(name)(paramSet);
    paramSet.reportUnused();

    return si;
  }

  VolumeIntegrator _makeVolumeIntegrator(String name, ParamSet paramSet) {
    if (overrides != null && overrides.volumeIntegratorName != null) {
      name = overrides.volumeIntegratorName;
      paramSet = overrides.volumeIntegratorParams;
    }

    if (Plugin.volumeIntegrator(name) == null) {
      LogWarning('Volume Integrator \'${name}\' unknown.');
      return null;
    }

    VolumeIntegrator vi = Plugin.volumeIntegrator(name)(paramSet);
    paramSet.reportUnused();

    return vi;
  }

  Primitive _makeAccelerator(String name, List<Primitive> prims,
                            ParamSet paramSet) {
    if (overrides != null && overrides.acceleratorName != null) {
      name = overrides.acceleratorName;
      paramSet = overrides.acceleratorParams;
    }

    if (Plugin.accelerator(name) == null) {
      LogWarning('Accelerator \'${name}\' unknown.');
      return null;
    }

    Primitive a = Plugin.accelerator(name)(prims, paramSet);
    paramSet.reportUnused();

    return a;
  }

  Camera _makeCamera(String name, ParamSet paramSet,
          TransformSet cam2worldSet, double transformStart,
          double transformEnd, Film film) {
    if (overrides != null && overrides.cameraName != null) {
      name = overrides.cameraName;
      paramSet = overrides.cameraParams;
    }

    if (Plugin.camera(name) == null) {
      LogWarning('Camera \'${name}\' unknown.');
      return null;
    }

    Transform cam2world0 = new Transform.from(cam2worldSet[0]);
    Transform cam2world1 = new Transform.from(cam2worldSet[1]);

    AnimatedTransform animatedCam2World = new AnimatedTransform(cam2world0,
        transformStart, cam2world1, transformEnd);

    Camera c = Plugin.camera(name)(paramSet, animatedCam2World, film);
    paramSet.reportUnused();

    return c;
  }

  PixelSampler _makePixelSampler(String name, ParamSet paramSet) {
    if (overrides != null && overrides.pixelSamplerName != null) {
      name = overrides.pixelSamplerName;
      paramSet = overrides.pixelSamplerParams;
    }

    if (Plugin.pixelSampler(name) == null) {
      LogWarning('PixelSampler \'${name}\' unknown.');
      return null;
    }

    PixelSampler s = Plugin.pixelSampler(name)(paramSet);
    paramSet.reportUnused();

    return s;
  }

  Sampler _makeSampler(String name, ParamSet paramSet, Film film,
                       Camera camera, PixelSampler pixels) {
    if (overrides != null && overrides.samplerName != null) {
      name = overrides.samplerName;
      paramSet = overrides.samplerParams;
    }

    if (Plugin.sampler(name) == null) {
      LogWarning('Sampler \'${name}\' unknown.');
      return null;
    }

    List<int> extent = [0, 0, 0, 0];
    film.getSampleExtent(extent);
    int w = extent[1] - extent[0];
    int h = extent[3] - extent[2];

    GetSubWindow(w, h, _renderOptions.taskNum, _renderOptions.taskCount,
                 extent);

    w = extent[1] - extent[0];
    h = extent[3] - extent[2];

    LogInfo('SAMPLER: $name [${extent[0]}, ${extent[2]}, $w, $h]');

    Sampler s = Plugin.sampler(name)(paramSet, extent[0], extent[2],
                                     w, h, camera, pixels);

    paramSet.reportUnused();

    return s;
  }

  Filter _makeFilter(String name, ParamSet paramSet) {
    if (overrides != null && overrides.filterName != null) {
      name = overrides.filterName;
      paramSet = overrides.filterParams;
    }

    if (Plugin.filter(name) == null) {
      LogWarning('Filter \'${name}\' unknown.');
      return null;
    }

    Filter f = Plugin.filter(name)(paramSet);
    paramSet.reportUnused();

    return f;
  }

  Film _makeFilm(String name, ParamSet paramSet, Filter filter,
                [PreviewCallback previewCallback]) {
    if (overrides != null && overrides.filmName != null) {
      name = overrides.filmName;
      paramSet = overrides.filmParams;
    }

    if (Plugin.film(name) == null) {
      LogWarning('Film \'${name}\' unknown.');
      return null;
    }

    Film f = Plugin.film(name)(paramSet, filter, previewCallback);

    paramSet.reportUnused();

    return f;
  }

  List _renders = [];
}
