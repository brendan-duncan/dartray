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
part of materials;

class MixMaterial extends Material {
  MixMaterial(this.m1, this.m2, this.scale);

  static MixMaterial Create(Transform xform,
          TextureParams mp, Material m1, Material m2) {
    Texture scale = mp.getSpectrumTexture("amount", new Spectrum(0.5));
    return new MixMaterial(m1, m2, scale);
  }

  BSDF getBSDF(DifferentialGeometry dgGeom,
               DifferentialGeometry dgShading) {
    BSDF b1 = m1.getBSDF(dgGeom, dgShading);
    BSDF b2 = m2.getBSDF(dgGeom, dgShading);
    Spectrum s1 = scale.evaluate(dgShading).clamp();
    Spectrum s2 = (new Spectrum(1.0) - s1).clamp();
    int n1 = b1.numComponents();
    int n2 = b2.numComponents();
    for (int i = 0; i < n1; ++i) {
      b1.bxdfs[i] = new ScaledBxDF(b1.bxdfs[i], s1);
    }
    for (int i = 0; i < n2; ++i) {
      b1.add(new ScaledBxDF(b2.bxdfs[i], s2));
    }
    return b1;
  }

  Material m1, m2;
  Texture scale;
}
