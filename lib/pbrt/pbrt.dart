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
library pbrt;

import 'package:image/image.dart';

import '../accelerators/accelerators.dart';
import '../core/core.dart';
import '../cameras/cameras.dart';
import '../film/film.dart';
import '../filters/filters.dart';
import '../lights/lights.dart';
import '../materials/materials.dart';
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
part 'render_options.dart';
part 'transform_set.dart';

typedef Aggregate AcceleratorCreator(List<Primitive> prims, ParamSet ps);

typedef Camera CameraCreator(ParamSet params, AnimatedTransform cam2world,
                             Film film);

typedef Film FilmCreator(ParamSet params, Filter filter,
                         [Image image, PreviewCallback previewCallback]);

typedef Filter FilterCreator(ParamSet ps);

typedef SurfaceIntegrator SurfaceIntegratorCreator(ParamSet ps);

typedef VolumeIntegrator VolumeIntegratorCreator(ParamSet ps);

typedef Light LightCreator(Transform light2world, ParamSet paramSet);

typedef Light AreaLightCreator(Transform light2world, ParamSet paramSet,
                               Shape shape);

typedef Material MaterialCreator(Transform xform, TextureParams mp);

typedef Sampler SamplerCreator(ParamSet params, Film film, Camera camera);

typedef Shape ShapeCreator(Transform o2w, Transform w2o,
                           bool reverseOrientation, ParamSet params);

typedef Texture TextureCreator(Transform tex2world, TextureParams tp);

typedef VolumeRegion VolumeRegionCreator(Transform volume2world,
                                         ParamSet params);



/**
 * Pbrt provides an API and file format compatible with the PBRT
 * rendering system API and format.
 */
class Pbrt {
  static Map<String, AcceleratorCreator> _accelerators = {};
  static Map<String, CameraCreator> _cameras = {};
  static Map<String, FilmCreator> _films = {};
  static Map<String, FilterCreator> _filters = {};
  static Map<String, SurfaceIntegratorCreator> _surfaceIntegrators = {};
  static Map<String, VolumeIntegratorCreator> _volumeIntegrators = {};
  static Map<String, LightCreator> _lights = {};
  static Map<String, AreaLightCreator> _areaLights = {};
  static Map<String, MaterialCreator> _materials = {};
  static Map<String, SamplerCreator> _samplers = {};
  static Map<String, ShapeCreator> _shapes = {};
  static Map<String, TextureCreator> _floatTextures = {};
  static Map<String, TextureCreator> _spectrumTextures = {};
  static Map<String, VolumeRegionCreator> _volumeRegions = {};

  static void registerAccelerator(String name, AcceleratorCreator func) {
    _accelerators[name] = func;
  }
  static void registerCamera(String name, CameraCreator func) {
    _cameras[name] = func;
  }
  static void registerFilm(String name, FilmCreator func) {
    _films[name] = func;
  }
  static void registerFilter(String name, FilterCreator func) {
    _filters[name] = func;
  }
  static void registerSurfaceIntegrator(String name,
                                        SurfaceIntegratorCreator func) {
    _surfaceIntegrators[name] = func;
  }
  static void registerVolumeIntegrator(String name,
                                        VolumeIntegratorCreator func) {
    _volumeIntegrators[name] = func;
  }
  static void registerLight(String name, LightCreator func) {
    _lights[name] = func;
  }
  static void registerAreaLight(String name, AreaLightCreator func) {
    _areaLights[name] = func;
  }
  static void registerMaterial(String name, MaterialCreator func) {
    _materials[name] = func;
  }
  static void registerSampler(String name, SamplerCreator func) {
    _samplers[name] = func;
  }
  static void registerShape(String name, ShapeCreator func) {
    _shapes[name] = func;
  }
  static void registerFloatTexture(String name, TextureCreator func) {
    _floatTextures[name] = func;
  }
  static void registerSpectrumTexture(String name, TextureCreator func) {
    _spectrumTextures[name] = func;
  }
  static void registerVolumeRegion(String name, VolumeRegionCreator func) {
    _volumeRegions[name] = func;
  }

  static void _registerStandardNodes() {
    if (_shapes.containsKey('sphere')) {
      return;
    }

    registerAccelerator('bvh', BVHAccel.Create);
    registerAccelerator('grid', GridAccel.Create);
    registerAccelerator('kdtree', KdTreeAccel.Create);
    registerAccelerator('naive', NaiveAccel.Create);

    registerCamera('environment', EnvironmentCamera.Create);
    registerCamera('orthographic', OrthographicCamera.Create);
    registerCamera('perspective', PerspectiveCamera.Create);

    registerFilm('image', ImageFilm.Create);

    registerFilter('box', BoxFilter.Create);
    registerFilter('gaussian', GaussianFilter.Create);
    registerFilter('sinc', LanczosSincFilter.Create);
    registerFilter('mitchell', MitchellFilter.Create);
    registerFilter('triangle', TriangleFilter.Create);

    registerSurfaceIntegrator('ambientocclusion', AmbientOcclusionIntegrator.Create);
    registerSurfaceIntegrator('diffuseprt', DiffusePRTIntegrator.Create);
    registerSurfaceIntegrator('directlighting', DirectLightingIntegrator.Create);
    registerSurfaceIntegrator('glossyprt', GlossyPRTIntegrator.Create);
    registerSurfaceIntegrator('igi', IGIIntegrator.Create);
    registerSurfaceIntegrator('irradiancecache', IrradianceCacheIntegrator.Create);
    registerSurfaceIntegrator('path', PathIntegrator.Create);
    registerSurfaceIntegrator('photonmap', PhotonMapIntegrator.Create);
    // alias for photonmap
    registerSurfaceIntegrator('exphotonmap', PhotonMapIntegrator.Create);
    registerSurfaceIntegrator('whitted', WhittedIntegrator.Create);

    registerLight('distant', DistantLight.Create);
    registerLight('point', PointLight.Create);
    registerLight('spot', SpotLight.Create);

    registerAreaLight('diffuse', DiffuseAreaLight.Create);
    // alias for diffuse
    registerAreaLight('area', DiffuseAreaLight.Create);

    registerMaterial('glass', GlassMaterial.Create);
    registerMaterial('kdsubsurface', KdSubsurfaceMaterial.Create);
    registerMaterial('matte', MatteMaterial.Create);
    registerMaterial('metal', MetalMaterial.Create);
    registerMaterial('mirror', MirrorMaterial.Create);
    registerMaterial('plastic', PlasticMaterial.Create);
    registerMaterial('shinymetal', ShinyMetalMaterial.Create);
    registerMaterial('substrate', SubstrateMaterial.Create);
    registerMaterial('subsurface', SubsurfaceMaterial.Create);
    registerMaterial('translucent', TranslucentMaterial.Create);
    registerMaterial('uber', UberMaterial.Create);

    registerSampler('adaptive', AdaptiveSampler.Create);
    registerSampler('bestcandidate', BestCandidateSampler.Create);
    registerSampler('halton', HaltonSampler.Create);
    registerSampler('lowdiscrepancy', LowDiscrepancySampler.Create);
    registerSampler('random', RandomSampler.Create);
    registerSampler('stratified', StratifiedSampler.Create);

    registerShape('cone', ConeShape.Create);
    registerShape('cylinder', CylinderShape.Create);
    registerShape('disk', DiskShape.Create);
    registerShape('hyperboloid', HyperboloidShape.Create);
    registerShape('loopsubdiv', LoopSubdivisionShape.Create);
    registerShape('paraboloid', ParaboloidShape.Create);
    registerShape('sphere', SphereShape.Create);
    registerShape('trianglemesh', TriangleMeshShape.Create);

    registerFloatTexture('bilerp', BilerpTexture.CreateFloat);
    registerSpectrumTexture('bilerp', BilerpTexture.CreateSpectrum);
    registerFloatTexture('checkerboard', CheckerboardTexture.CreateFloat);
    registerSpectrumTexture('checkerboard', CheckerboardTexture.CreateSpectrum);
    registerFloatTexture('dots', DotsTexture.CreateFloat);
    registerSpectrumTexture('dots', DotsTexture.CreateSpectrum);
    registerFloatTexture('fbm', FBmTexture.CreateFloat);
    registerSpectrumTexture('fbm', FBmTexture.CreateSpectrum);
    registerFloatTexture('marble', MarbleTexture.CreateFloat);
    registerSpectrumTexture('marble', MarbleTexture.CreateSpectrum);
    registerFloatTexture('mix', MixTexture.CreateFloat);
    registerSpectrumTexture('mix', MixTexture.CreateSpectrum);
    registerFloatTexture('scale', ScaleTexture.CreateFloat);
    registerSpectrumTexture('scale', ScaleTexture.CreateSpectrum);
    registerFloatTexture('uv', UVTexture.CreateFloat);
    registerSpectrumTexture('uv', UVTexture.CreateSpectrum);
    registerFloatTexture('windy', WindyTexture.CreateFloat);
    registerSpectrumTexture('windy', WindyTexture.CreateSpectrum);
    registerFloatTexture('wrinkled', WrinkledTexture.CreateFloat);
    registerSpectrumTexture('wrinkled', WrinkledTexture.CreateSpectrum);

    registerVolumeIntegrator('emission', EmissionIntegrator.Create);

    registerVolumeRegion('homogenous', HomogeneousVolumeDensityRegion.Create);
  }

  Pbrt() {
    _registerStandardNodes();
  }

  OutputImage renderScene(String scene, [Image output]) {
    if (output != null) {
      setOutputImage(output);
    }

    if (!loadScene(scene)) {
      return null;
    }

    return _renderer.render(_scene);
  }

  Renderer getRenderer() => _renderer;

  Scene getScene() => _scene;

  bool loadScene(String file) {
    PbrtParser parser = new PbrtParser(this);
    parser.parse(file);
    return _renderer != null && _scene != null;
  }

  OutputImage render() {
    if (_renderer == null || _scene == null) {
      return null;
    }
    return _renderer.render(_scene);
  }

  static const int _MAX_TRANSFORMS = 2;
  static const int _START_TRANSFORM_BITS = (1 << 0);
  static const int _END_TRANSFORM_BITS = (1 << 1);
  static const int _ALL_TRANSFORMS_BITS = ((1 << _MAX_TRANSFORMS) - 1);

  static const int STATE_UNINITIALIZED = 0;
  static const int STATE_OPTIONS_BLOCK = 1;
  static const int STATE_WORLD_BLOCK = 2;

  int _currentApiState = STATE_OPTIONS_BLOCK;
  TransformSet _curTransform = new TransformSet();
  int _activeTransformBits = _ALL_TRANSFORMS_BITS;
  Map<String, TransformSet> _namedCoordinateSystems = {};
  RenderOptions _renderOptions = new RenderOptions();
  GraphicsState _graphicsState = new GraphicsState();
  List<GraphicsState> _pushedGraphicsStates = [];
  List<TransformSet> _pushedTransforms = [];
  List<int> _pushedActiveTransformBits = [];


  // API Initialization
  void init() {
  }

  // API Cleanup
  void cleanup() {
    //ProbesCleanup();
    if (_currentApiState == STATE_UNINITIALIZED) {
      LogWarning("cleanup() called without init().");
    } else if (_currentApiState == STATE_WORLD_BLOCK) {
      LogWarning("cleanup() called while inside world block.");
    }

    _currentApiState = STATE_UNINITIALIZED;
  }

  void setTask(int taskNum, int taskCount) {
    _renderOptions.taskNum = taskNum;
    _renderOptions.taskCount = taskCount;
  }

  void identity() {
    for (int i = 0; i < _MAX_TRANSFORMS; ++i) {
      if (_activeTransformBits & (1 << i) != 0) {
        _curTransform[i] = new Transform();
      }
    }
  }

  void translate(double dx, double dy, double dz) {
    for (int i = 0; i < _MAX_TRANSFORMS; ++i) {
      if (_activeTransformBits & (1 << i) != 0) {
        _curTransform[i] = _curTransform[i] *
                          Transform.Translate(new Vector(dx, dy, dz));
      }
    }
  }

  void transform(Matrix4x4 tr) {
    for (int i = 0; i < _MAX_TRANSFORMS; ++i) {
      if (_activeTransformBits & (1 << i) != 0) {
        _curTransform[i] = new Transform(tr);
      }
    }
  }

  void concatTransform(Matrix4x4 tr) {
    for (int i = 0; i < _MAX_TRANSFORMS; ++i) {
      if (_activeTransformBits & (1 << i) != 0) {
        _curTransform[i] = _curTransform[i] * new Transform(tr);
      }
    }
  }

  void rotate(double angle, double dx, double dy, double dz) {
    for (int i = 0; i < _MAX_TRANSFORMS; ++i) {
      if (_activeTransformBits & (1 << i) != 0) {
        _curTransform[i] = _curTransform[i] *
                          Transform.Rotate(angle, new Vector(dx, dy, dz));
      }
    }
  }

  void scale(double sx, double sy, double sz) {
    for (int i = 0; i < _MAX_TRANSFORMS; ++i) {
      if (_activeTransformBits & (1 << i) != 0) {
        _curTransform[i] = _curTransform[i] * Transform.Scale(sx, sy, sz);
      }
    }
  }

  void lookAt(double ex, double ey, double ez, double lx, double ly,
              double lz, double ux, double uy, double uz) {
    for (int i = 0; i < _MAX_TRANSFORMS; ++i) {
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
      _curTransform = _namedCoordinateSystems[name];
    } else {
      LogWarning("Couldn't find named coordinate system \"$name\"");
    }
  }

  void activeTransformAll() {
    _activeTransformBits = _ALL_TRANSFORMS_BITS;
  }

  void activeTransformEndTime() {
    _activeTransformBits = _END_TRANSFORM_BITS;
  }

  void activeTransformStartTime() {
    _activeTransformBits = _START_TRANSFORM_BITS;
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

  void setOutputImage(Image output) {
    _renderOptions.outputImage = output;
  }

  void setPreviewCallback(PreviewCallback cb) {
    _renderOptions.previewCallback = cb;
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
    _renderOptions.surfIntegratorName = name;
    _renderOptions.surfIntegratorParams = params;
  }

  void volumeIntegrator(String name, ParamSet params) {
    _renderOptions.volIntegratorName = name;
    _renderOptions.volIntegratorParams = params;
  }

  void renderer(String name, ParamSet params) {
    _renderOptions.rendererName = name;
    _renderOptions.rendererParams = params;
  }

  void camera(String name, ParamSet params) {
    _renderOptions.cameraName = name;
    _renderOptions.cameraParams = params;
    _renderOptions.cameraToWorld = TransformSet.Inverse(_curTransform);
    _namedCoordinateSystems["camera"] = _renderOptions.cameraToWorld;
  }

  void worldBegin() {
    _currentApiState = STATE_WORLD_BLOCK;
    for (int i = 0; i < _MAX_TRANSFORMS; ++i) {
      _curTransform[i] = new Transform();
    }
    _activeTransformBits = _ALL_TRANSFORMS_BITS;
    _namedCoordinateSystems["world"] = _curTransform;
  }

  void attributeBegin() {
    _pushedGraphicsStates.add(new GraphicsState.from(_graphicsState));
    _pushedTransforms.add(new TransformSet.from(_curTransform));
    _pushedActiveTransformBits.add(_activeTransformBits);
  }

  void attributeEnd() {
    if (_pushedGraphicsStates.isEmpty) {
      LogWarning("Unmatched attributeEnd() encountered. "
              "Ignoring it.");
      return;
    }
    _graphicsState = _pushedGraphicsStates.removeLast();
    _curTransform = _pushedTransforms.removeLast();
  }

  void transformBegin() {
    _pushedTransforms.add(_curTransform);
    _pushedActiveTransformBits.add(_activeTransformBits);
  }

  void transformEnd() {
    if (_pushedTransforms.isNotEmpty) {
      LogWarning("Unmatched pbrtTransformEnd() encountered. "
          "Ignoring it.");
      return;
    }

    _curTransform = _pushedTransforms.removeLast();
    _activeTransformBits = _pushedActiveTransformBits.removeLast();
  }

  void texture(String name, String type, String texname, ParamSet params) {
    TextureParams tp = new TextureParams(params, params,
                            _graphicsState.doubleTextures,
                            _graphicsState.spectrumTextures);
    if (type == "float")  {
      // Create _double_ texture and store in _doubleTextures_
      if (_graphicsState.doubleTextures.containsKey(name)) {
        LogWarning("Texture \"$name\" being redefined");
      }

      //WARN_IF_ANIMATED_TRANSFORM("Texture");
      Texture ft = _makeFloatTexture(texname, _curTransform[0], tp);
      if (ft != null) {
        _graphicsState.doubleTextures[name] = ft;
      }
    } else if (type == "color" || type == "spectrum")  {
      // Create _color_ texture and store in _spectrumTextures_
      if (_graphicsState.spectrumTextures.containsKey(name)) {
        LogWarning("Texture \"$name\" being redefined");
      }

      //WARN_IF_ANIMATED_TRANSFORM("Texture");
      Texture st = _makeSpectrumTexture(texname, _curTransform[0], tp);
      if (st != null) {
        _graphicsState.spectrumTextures[name] = st;
      }
    } else {
      LogWarning("Texture type \"$type\" unknown.");
    }
  }

  void material(String name, ParamSet params) {
    _graphicsState.material = name;
    _graphicsState.materialParams = params;
    _graphicsState.currentNamedMaterial = "";
  }

  void makeNamedMaterial(String name, ParamSet params) {
    // error checking, warning if replace, what to use for transform?
    TextureParams mp = new TextureParams(params, _graphicsState.materialParams,
                     _graphicsState.doubleTextures,
                     _graphicsState.spectrumTextures);
    String matName = mp.findString("type");

    //WARN_IF_ANIMATED_TRANSFORM("MakeNamedMaterial");
    if (matName.isEmpty) {
      LogWarning("No parameter string \"type\" found in MakeNamedMaterial");
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
    //WARN_IF_ANIMATED_TRANSFORM("LightSource");
    Light lt = _makeLight(name, _curTransform[0], params);
    if (lt == null) {
      LogWarning("pbrtLightSource: light type \"$name\" unknown.");
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
      if (_graphicsState.areaLight != "") {
        area = _makeAreaLight(_graphicsState.areaLight, _curTransform[0],
                             _graphicsState.areaLightParams, shape);
      }

      prim = new GeometricPrimitive(shape, mtl, area);
    } else {
      // Create primitive for animated shape

      // Create initial _Shape_ for animated shape
      if (_graphicsState.areaLight != "") {
        LogWarning("Ignoring currently set area light when creating "
              "animated shape");
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
      assert(_MAX_TRANSFORMS == 2);

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
        LogWarning("Area lights not supported with object instancing");
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
      LogWarning("ObjectBegin called inside of instance definition");
    }
    _renderOptions.instances[name] = new List<Primitive>();
    _renderOptions.currentInstance = _renderOptions.instances[name];
  }


  void objectEnd() {
    if (_renderOptions.currentInstance == null) {
      LogWarning("ObjectEnd called outside of instance definition");
    }
    _renderOptions.currentInstance = null;
    attributeEnd();
  }


  void objectInstance(String name) {
    // Object instance error checking
    if (_renderOptions.currentInstance != null) {
      LogWarning("ObjectInstance can't be called inside instance definition");
      return;
    }

    if (!_renderOptions.instances.containsKey(name)) {
      LogWarning("Unable to find instance named \"$name\"");
      return;
    }

    List<Primitive> inst = _renderOptions.instances[name];

    if (inst.length == 0) {
      return;
    }
    if (inst.length > 1 || !inst[0].canIntersect()) {
      // Refine instance _Primitive_s and create aggregate
      Primitive accel = _makeAccelerator(_renderOptions.acceleratorName,
                           inst, _renderOptions.acceleratorParams);

      if (accel == null) {
        accel = _makeAccelerator("bvh", inst, new ParamSet());
      }

      if (accel == null) {
        LogSevere("Unable to create \"bvh\" accelerator");
      }

      inst.clear();
      inst.add(accel);
    }

    assert(_MAX_TRANSFORMS == 2);

    Transform world2instance0 = Transform.Inverse(_curTransform[0]);
    Transform world2instance1 = Transform.Inverse(_curTransform[1]);

    AnimatedTransform animatedWorldToInstance = new AnimatedTransform(
        world2instance0, _renderOptions.transformStartTime,
        world2instance1, _renderOptions.transformEndTime);

    Primitive prim = new TransformedPrimitive(inst[0], animatedWorldToInstance);
    _renderOptions.primitives.add(prim);
  }


  void worldEnd() {
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
    _renderer = _makeRenderer();
    _scene = _makeScene();

    // Clean up after rendering
    _graphicsState = new GraphicsState();

    _currentApiState = STATE_OPTIONS_BLOCK;
    //ProbesPrint(stdout);
    for (int i = 0; i < _MAX_TRANSFORMS; ++i) {
      _curTransform.t[i] = new Transform();
    }

    _activeTransformBits = _ALL_TRANSFORMS_BITS;
    _namedCoordinateSystems.clear();
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
        _renderOptions.primitives, _renderOptions.acceleratorParams);
    if (accelerator == null) {
      accelerator = _makeAccelerator('bvh', _renderOptions.primitives,
                                     new ParamSet());
    }

    if (accelerator == null) {
      LogSevere("Unable to create \"bvh\" accelerator.");
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

    Film film = _makeFilm(_renderOptions.filmName, _renderOptions.filmParams,
                          filter, _renderOptions.outputImage,
                          _renderOptions.previewCallback);
    if (film == null) {
      LogSevere("Unable to create film.");
    }

    Camera camera = _makeCamera(_renderOptions.cameraName,
                                _renderOptions.cameraParams,
                                _renderOptions.cameraToWorld,
                                _renderOptions.transformStartTime,
                                _renderOptions.transformEndTime,
                                film);

    if (camera == null) {
      LogSevere("Unable to create camera.");
    }

    /*if (RendererName == "metropolis") {
      renderer = MetropolisRenderer.Create(RendererParams, camera);
      RendererParams.reportUnused();
      // Warn if no light sources are defined
      if (lights.length == 0) {
        Warning("No light sources defined in scene; "
                "possibly rendering a black image.");
      }
    } else if (RendererName == "createprobes") {
      // Create surface and volume integrators
      SurfaceIntegrator *surfaceIntegrator = MakeSurfaceIntegrator(SurfIntegratorName,
          SurfIntegratorParams);
      if (!surfaceIntegrator) Severe("Unable to create surface integrator.");
      VolumeIntegrator *volumeIntegrator = MakeVolumeIntegrator(VolIntegratorName,
          VolIntegratorParams);
      if (!volumeIntegrator) Severe("Unable to create volume integrator.");
      renderer = CreateRadianceProbesRenderer(camera, surfaceIntegrator, volumeIntegrator, RendererParams);
      RendererParams.ReportUnused();
      // Warn if no light sources are defined
      if (lights.size() == 0)
          Warning("No light sources defined in scene; "
              "possibly rendering a black image.");
    } else if (RendererName == "aggregatetest") {
      renderer = CreateAggregateTestRenderer(RendererParams, primitives);
      RendererParams.ReportUnused();
    } else if (RendererName == "surfacepoints") {
      Point pCamera = camera.CameraToWorld(camera.shutterOpen, Point(0, 0, 0));
      renderer = CreateSurfacePointsRenderer(RendererParams, pCamera, camera.shutterOpen);
      RendererParams.ReportUnused();
    } else*/ {
      if (_renderOptions.rendererName != "sampler") {
        LogWarning("Renderer type \"${_renderOptions.rendererName}\" unknown.  Using \"sampler\".");
      }

      _renderOptions.rendererParams.reportUnused();

      Sampler sampler = _makeSampler(_renderOptions.samplerName,
                                     _renderOptions.samplerParams,
                                     camera.film, camera);
      if (sampler == null) {
        LogSevere("Unable to create sampler.");
      }

      // Create surface and volume integrators
      SurfaceIntegrator surfaceIntegrator =
          _makeSurfaceIntegrator(_renderOptions.surfIntegratorName,
                                 _renderOptions.surfIntegratorParams);
      if (surfaceIntegrator == null) {
        LogSevere("Unable to create surface integrator.");
      }

      VolumeIntegrator volumeIntegrator =
          _makeVolumeIntegrator(_renderOptions.volIntegratorName,
                                _renderOptions.volIntegratorParams);
      if (volumeIntegrator == null) {
        LogSevere("Unable to create volume integrator.");
      }

      renderer = new SamplerRenderer(sampler, camera, surfaceIntegrator,
                                     volumeIntegrator, _renderOptions.taskNum,
                                     _renderOptions.taskCount);

      // Warn if no light sources are defined
      if (_renderOptions.lights.length == 0) {
        LogWarning("No light sources defined in scene; "
              "possibly rendering a black image.");
      }
    }

    return renderer;
  }

  Shape _makeShape(String name, Transform object2world, Transform world2object,
                  bool reverseOrientation, ParamSet paramSet) {
    if (!_shapes.containsKey(name)) {
      LogWarning('Shape \'${name}\' unknown.');
      return null;
    }

    Shape s = _shapes[name](object2world, world2object, reverseOrientation,
                            paramSet);
    paramSet.reportUnused();

    return s;
  }

  Material _createMaterial(ParamSet params) {
    TextureParams mp = new TextureParams(params, _graphicsState.materialParams,
        _graphicsState.doubleTextures,
        _graphicsState.spectrumTextures);

    Material mtl;
    if (_graphicsState.currentNamedMaterial != "" &&
        _graphicsState.namedMaterials.containsKey(_graphicsState.currentNamedMaterial)) {
      mtl = _graphicsState.namedMaterials[_graphicsState.currentNamedMaterial];
    }

    if (mtl == null) {
      mtl = _makeMaterial(_graphicsState.material, _curTransform[0], mp);
    }

    if (mtl == null) {
      mtl = _makeMaterial("matte", _curTransform[0], mp);
    }

    if (mtl == null) {
      LogSevere("Unable to create \"matte\" material?!");
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

    if (!_materials.containsKey(name)) {
      LogWarning('Material \'$name\' unknown.');
      return null;
    }

    Material material = _materials[name](mtl2world, mp);
    mp.reportUnused();

    return material;
  }

  Texture _makeFloatTexture(String name, Transform tex2world, TextureParams tp) {
    if (!_floatTextures.containsKey(name)) {
      LogWarning('Texture \'${name}\' unknown.');
      return ConstantTexture.CreateFloat(tex2world, tp);
    }

    Texture t = _floatTextures[name](tex2world, tp);
    tp.reportUnused();

    return t;
  }

  Texture _makeSpectrumTexture(String name, Transform tex2world,
                              TextureParams tp) {
    if (!_spectrumTextures.containsKey(name)) {
      LogWarning('Texture \'${name}\' unknown.');
      return ConstantTexture.CreateSpectrum(tex2world, tp);
    }

    Texture t = _spectrumTextures[name](tex2world, tp);
    tp.reportUnused();

    return t;
  }

  Light _makeLight(String name, Transform light2world, ParamSet paramSet) {
    if (!_lights.containsKey(name)) {
      LogWarning('Light \'${name}\' unknown.');
      return null;
    }

    Light l = _lights[name](light2world, paramSet);
    paramSet.reportUnused();

    return l;
  }

  AreaLight _makeAreaLight(String name, Transform light2world, ParamSet paramSet,
                          Shape shape) {
    if (!_areaLights.containsKey(name)) {
      LogWarning('Area Light \'${name}\' unknown.');
      return null;
    }

    Light l = _areaLights[name](light2world, paramSet, shape);
    paramSet.reportUnused();

    return l;
  }

  VolumeRegion _makeVolumeRegion(String name, Transform volume2world,
                                ParamSet paramSet) {
    if (!_volumeRegions.containsKey(name)) {
      LogWarning('Volume Region \'${name}\' unknown.');
      return null;
    }

    VolumeRegion v = _volumeRegions[name](volume2world, paramSet);
    paramSet.reportUnused();

    return v;
  }

  SurfaceIntegrator _makeSurfaceIntegrator(String name, ParamSet paramSet) {
    if (!_surfaceIntegrators.containsKey(name)) {
      LogWarning('Surface Integrator \'${name}\' unknown.');
      return null;
    }

    SurfaceIntegrator si = _surfaceIntegrators[name](paramSet);
    paramSet.reportUnused();

    return si;
  }

  VolumeIntegrator _makeVolumeIntegrator(String name, ParamSet paramSet) {
    if (!_volumeIntegrators.containsKey(name)) {
      LogWarning('Volume Integrator \'${name}\' unknown.');
      return null;
    }

    VolumeIntegrator vi = _volumeIntegrators[name](paramSet);
    paramSet.reportUnused();

    return vi;
  }

  Primitive _makeAccelerator(String name, List<Primitive> prims,
                            ParamSet paramSet) {
    if (!_accelerators.containsKey(name)) {
      LogWarning('Accelerator \'${name}\' unknown.');
      return null;
    }

    Primitive a = _accelerators[name](prims, paramSet);
    paramSet.reportUnused();

    return a;
  }

  Camera _makeCamera(String name, ParamSet paramSet,
          TransformSet cam2worldSet, double transformStart,
          double transformEnd, Film film) {
    if (!_cameras.containsKey(name)) {
      LogWarning('Camera \'${name}\' unknown.');
      return null;
    }

    Transform cam2world0 = new Transform.from(cam2worldSet[0]);
    Transform cam2world1 = new Transform.from(cam2worldSet[1]);

    AnimatedTransform animatedCam2World = new AnimatedTransform(cam2world0,
        transformStart, cam2world1, transformEnd);

    Camera c = _cameras[name](paramSet, animatedCam2World, film);
    paramSet.reportUnused();

    return c;
  }

  Sampler _makeSampler(String name,
          ParamSet paramSet, Film film, Camera camera) {
    if (!_samplers.containsKey(name)) {
      LogWarning('Sampler \'${name}\' unknown.');
      return null;
    }

    Sampler s = _samplers[name](paramSet, film, camera);
    paramSet.reportUnused();

    return s;
  }

  Filter _makeFilter(String name, ParamSet paramSet) {
    if (!_filters.containsKey(name)) {
      LogWarning('Filter \'${name}\' unknown.');
      return null;
    }

    Filter f = _filters[name](paramSet);
    paramSet.reportUnused();

    return f;
  }

  Film _makeFilm(String name, ParamSet paramSet, Filter filter,
                [Image outputImage, PreviewCallback previewCallback]) {
    if (!_films.containsKey(name)) {
      LogWarning('Film \'${name}\' unknown.');
      return null;
    }

    Film f = _films[name](paramSet, filter, outputImage, previewCallback);

    paramSet.reportUnused();

    return f;
  }

  Renderer _renderer;
  Scene _scene;
}
