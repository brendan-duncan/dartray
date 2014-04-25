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

class GraphicsState {
  GraphicsState() {
    material = 'matte';
    reverseOrientation = false;
  }

  GraphicsState.from(GraphicsState other) :
    doubleTextures = new Map.from(other.doubleTextures),
    spectrumTextures = new Map.from(other.spectrumTextures),
    materialParams = new ParamSet.from(other.materialParams),
    material = other.material,
    namedMaterials = new Map.from(other.namedMaterials),
    currentNamedMaterial = other.currentNamedMaterial,
    areaLightParams = new ParamSet.from(other.areaLightParams),
    areaLight = other.areaLight,
    reverseOrientation = other.reverseOrientation;

  Map<String, Texture> doubleTextures = {};
  Map<String, Texture> spectrumTextures = {};
  ParamSet materialParams = new ParamSet();
  String material;
  Map<String, Material> namedMaterials = {};
  String currentNamedMaterial;
  ParamSet areaLightParams = new ParamSet();
  String areaLight = '';
  bool reverseOrientation;
}
