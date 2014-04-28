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

typedef Aggregate AcceleratorCreator(List<Primitive> prims, ParamSet ps);

typedef Camera CameraCreator(ParamSet params, AnimatedTransform cam2world,
                             Film film);

typedef Film FilmCreator(ParamSet params, Filter filter,
                         [PreviewCallback previewCallback]);

typedef Filter FilterCreator(ParamSet ps);

typedef SurfaceIntegrator SurfaceIntegratorCreator(ParamSet ps);

typedef VolumeIntegrator VolumeIntegratorCreator(ParamSet ps);

typedef Light LightCreator(Transform light2world, ParamSet paramSet);

typedef Light AreaLightCreator(Transform light2world, ParamSet paramSet,
                               Shape shape);

typedef Material MaterialCreator(Transform xform, TextureParams mp);

typedef PixelSampler PixelSamplerCreator(ParamSet params);

typedef Sampler SamplerCreator(ParamSet params, int x, int y,
                               int width, int height, Camera camera,
                               PixelSampler pixels);

typedef Shape ShapeCreator(Transform o2w, Transform w2o,
                           bool reverseOrientation, ParamSet params);

typedef Texture TextureCreator(Transform tex2world, TextureParams tp);

typedef VolumeRegion VolumeRegionCreator(Transform volume2world,
                                         ParamSet params);

typedef Renderer RendererCreator(ParamSet params);



class Plugin {
  static AcceleratorCreator accelerator(String name) => _accelerators[name];

  static CameraCreator camera(String name) => _cameras[name];

  static FilmCreator film(String name) => _films[name];

  static FilterCreator filter(String name) => _filters[name];

  static SurfaceIntegratorCreator surfaceIntegrator(String name) =>
      _surfaceIntegrators[name];

  static VolumeIntegratorCreator volumeIntegrator(String name) =>
      _volumeIntegrators[name];

  static LightCreator light(String name) => _lights[name];

  static AreaLightCreator areaLight(String name) => _areaLights[name];

  static MaterialCreator material(String name) => _materials[name];

  static SamplerCreator sampler(String name) => _samplers[name];

  static PixelSamplerCreator pixelSampler(String name) => _pixelSamplers[name];

  static ShapeCreator shape(String name) => _shapes[name];

  static TextureCreator floatTexture(String name) => _floatTextures[name];

  static TextureCreator spectrumTexture(String name) =>
      _spectrumTextures[name];

  static VolumeRegionCreator volumeRegion(String name) => _volumeRegions[name];

  static RendererCreator renderer(String name) => _renderers[name];

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

  static void registerPixelSampler(String name, PixelSamplerCreator func) {
    _pixelSamplers[name] = func;
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

  static void registerRenderer(String name, RendererCreator func) {
    _renderers[name] = func;
  }

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
  static Map<String, PixelSamplerCreator> _pixelSamplers = {};
  static Map<String, ShapeCreator> _shapes = {};
  static Map<String, TextureCreator> _floatTextures = {};
  static Map<String, TextureCreator> _spectrumTextures = {};
  static Map<String, VolumeRegionCreator> _volumeRegions = {};
  static Map<String, RendererCreator> _renderers = {};
}
