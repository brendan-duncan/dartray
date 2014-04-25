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

/**
 * Allows the client application to override the rendering options defined in
 * scene files.
 */
class RenderOverrides {
  bool quickRender = false;
  double resolutionScale = 1.0;
  int samplingMode = Sampler.FULL_SAMPLING;
  String filterName;
  ParamSet filterParams;
  String filmName;
  ParamSet filmParams;
  String pixelSamplerName;
  ParamSet pixelSamplerParams;
  String samplerName;
  ParamSet samplerParams;
  String acceleratorName;
  ParamSet acceleratorParams;
  String rendererName;
  ParamSet rendererParams;
  String surfaceIntegratorName;
  ParamSet surfaceIntegratorParams;
  String volumeIntegratorName;
  ParamSet volumeIntegratorParams;
  String cameraName;
  ParamSet cameraParams;

  static RenderOverrides global;

  static bool HasOverrides() {
    return global != null;
  }

  static bool QuickRender() {
    if (global != null) {
      return global.quickRender;
    }
    return false;
  }

  static double ResolutionScale() {
    if (global != null) {
      return global.resolutionScale;
    }
    return 1.0;
  }

  static int SamplingMode() {
    if (global != null) {
      return global.samplingMode;
    }
    return Sampler.FULL_SAMPLING;
  }

  RenderOverrides() {
    global = this;
  }

  RenderOverrides.fromJson(Map json) {
    global = this;
    if (json.containsKey('quickRender')) {
      quickRender = json['quickRender'];
    }
    if (json.containsKey('resolutionScale')) {
      resolutionScale = json['resolutionScale'];
    }
    if (json.containsKey('samplingMode')) {
      samplingMode = json['samplingMode'];
    }
    if (json.containsKey('filter')) {
      filterName = json['filter']['name'];
      filterParams = new ParamSet.fromJson(json['filter']['params']);
    }
    if (json.containsKey('film')) {
      filmName = json['film']['name'];
      filmParams = new ParamSet.fromJson(json['film']['params']);
    }
    if (json.containsKey('pixelSampler')) {
      pixelSamplerName = json['pixelSampler']['name'];
      pixelSamplerParams = new ParamSet.fromJson(json['pixelSampler']['params']);
    }
    if (json.containsKey('sampler')) {
      samplerName = json['sampler']['name'];
      samplerParams = new ParamSet.fromJson(json['sampler']['params']);
    }
    if (json.containsKey('accelerator')) {
      acceleratorName = json['accelerator']['name'];
      acceleratorParams = new ParamSet.fromJson(json['accelerator']['params']);
    }
    if (json.containsKey('renderer')) {
      rendererName = json['renderer']['name'];
      rendererParams = new ParamSet.fromJson(json['renderer']['params']);
    }
    if (json.containsKey('surfaceIntegrator')) {
      surfaceIntegratorName = json['surfaceIntegrator']['name'];
      surfaceIntegratorParams =
          new ParamSet.fromJson(json['surfaceIntegrator']['params']);
    }
    if (json.containsKey('volumeIntegrator')) {
      volumeIntegratorName = json['volumeIntegrator']['name'];
      volumeIntegratorParams =
          new ParamSet.fromJson(json['volumeIntegrator']['params']);
    }
    if (json.containsKey('camera')) {
      cameraName = json['camera']['name'];
      cameraParams = new ParamSet.fromJson(json['camera']['params']);
    }
  }

  Map toJson() {
    Map m = {};
    if (quickRender) {
      m['quickRender'] = quickRender;
    }
    if (resolutionScale != null) {
      m['resolutionScale'] = resolutionScale;
    }
    if (samplingMode != null) {
      m['samplingMode'] = samplingMode;
    }
    if (filterName != null) {
      m['filter'] = {'name': filterName, 'params': filterParams.toJson()};
    }
    if (filmName != null) {
      m['film'] = {'name': filmName, 'params': filmParams.toJson()};
    }
    if (pixelSamplerName != null) {
      m['pixelSampler'] = {'name': pixelSamplerName,
                           'params': pixelSamplerParams.toJson()};
    }
    if (samplerName != null) {
      m['sampler'] = {'name': samplerName, 'params': samplerParams.toJson()};
    }
    if (acceleratorName != null) {
      m['accelerator'] = {'name': acceleratorName,
                          'params': acceleratorParams.toJson()};
    }
    if (rendererName != null) {
      m['renderer'] = {'name': rendererName,
                       'params': rendererParams.toJson()};
    }
    if (surfaceIntegratorName != null) {
      m['surfaceIntegrator'] = {'name': surfaceIntegratorName,
                                'params': surfaceIntegratorParams.toJson()};
    }
    if (volumeIntegratorName != null) {
      m['volumeIntegrator'] = {'name': volumeIntegratorName,
                               'params': volumeIntegratorParams.toJson()};
    }
    if (cameraName != null) {
      m['camera'] = {'name': cameraName,
                     'params': cameraParams.toJson()};
    }
    return m;
  }

  void setFilter(String name, [params]) {
    if (params is Map) {
      params = new ParamSet.fromJson(params);
    } else if (params == null) {
      params = new ParamSet();
    }
    if (params is! ParamSet) {
      return;
    }
    filterName = name;
    filterParams = params;
  }

  void setFilm(String name, [params]) {
    if (params is Map) {
      params = new ParamSet.fromJson(params);
    } else if (params == null) {
      params = new ParamSet();
    }
    if (params is! ParamSet) {
      return;
    }
    filmName = name;
    filmParams = params;
  }

  void setPixelSampler(String name, [params]) {
    if (params is Map) {
      params = new ParamSet.fromJson(params);
    } else if (params == null) {
      params = new ParamSet();
    }
    if (params is! ParamSet) {
      return;
    }
    pixelSamplerName = name;
    pixelSamplerParams = params;
  }

  void setSampler(String name, [params]) {
    if (params is Map) {
      params = new ParamSet.fromJson(params);
    } else if (params == null) {
      params = new ParamSet();
    }
    if (params is! ParamSet) {
      return;
    }
    samplerName = name;
    samplerParams = params;
  }

  void setAccelerator(String name, [params]) {
    if (params is Map) {
      params = new ParamSet.fromJson(params);
    } else if (params == null) {
      params = new ParamSet();
    }
    if (params is! ParamSet) {
      return;
    }
    acceleratorName = name;
    acceleratorParams = params;
  }

  void setRenderer(String name, [params]) {
    if (params is Map) {
      params = new ParamSet.fromJson(params);
    } else if (params == null) {
      params = new ParamSet();
    }
    if (params is! ParamSet) {
      return;
    }
    rendererName = name;
    rendererParams = params;
  }

  void setSurfaceIntegrator(String name, [params]) {
    if (params is Map) {
      params = new ParamSet.fromJson(params);
    } else if (params == null) {
      params = new ParamSet();
    }
    if (params is! ParamSet) {
      return;
    }
    surfaceIntegratorName = name;
    surfaceIntegratorParams = params;
  }

  void setVolumeIntegrator(String name, [params]) {
    if (params is Map) {
      params = new ParamSet.fromJson(params);
    } else if (params == null) {
      params = new ParamSet();
    }
    if (params is! ParamSet) {
      return;
    }
    volumeIntegratorName = name;
    volumeIntegratorParams = params;
  }

  void setCamera(String name, [params]) {
    if (params is Map) {
      params = new ParamSet.fromJson(params);
    } else if (params == null) {
      params = new ParamSet();
    }
    if (params is! ParamSet) {
      return;
    }
    cameraName = name;
    cameraParams = params;
  }
}
