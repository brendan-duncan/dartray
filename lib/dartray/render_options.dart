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

class RenderOptions {
  RenderOptions() {
    transformStartTime = 0.0;
    transformEndTime = 1.0;
    filterName = 'box';
    filmName = 'image';
    pixelSamplerName = 'tile';
    samplerName = 'lowdiscrepancy';
    acceleratorName = 'bvh';
    rendererName = 'sampler';
    surfaceIntegratorName = 'directlighting';
    volumeIntegratorName = 'emission';
    cameraName = 'perspective';
    currentInstance = null;
    taskNum = 0;
    taskCount = 1;
  }

  double transformStartTime;
  double transformEndTime;
  String filterName;
  ParamSet filterParams = new ParamSet();
  String filmName;
  ParamSet filmParams = new ParamSet();
  PreviewCallback previewCallback;
  String pixelSamplerName;
  ParamSet pixelSamplerParams = new ParamSet();
  String samplerName;
  ParamSet samplerParams = new ParamSet();
  String acceleratorName;
  ParamSet acceleratorParams = new ParamSet();
  String rendererName;
  ParamSet rendererParams = new ParamSet();
  String surfaceIntegratorName;
  ParamSet surfaceIntegratorParams = new ParamSet();
  String volumeIntegratorName;
  ParamSet volumeIntegratorParams = new ParamSet();
  String cameraName;
  ParamSet cameraParams = new ParamSet();
  TransformSet cameraToWorld;
  List<Light> lights = [];
  List<Primitive> primitives = [];
  List<VolumeRegion> volumeRegions = [];
  Map<String, List<Primitive>> instances = {};
  List<Primitive> currentInstance = [];
  int taskNum;
  int taskCount;
}
