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

class GlassMaterial extends Material {
  GlassMaterial(this.Kr, this.Kt, this.index, this.bumpMap);

  static GlassMaterial Create(Transform xform, TextureParams mp) {
    Texture Kr = mp.getSpectrumTexture("Kr", new Spectrum(1.0));
    Texture Kt = mp.getSpectrumTexture("Kt", new Spectrum(1.0));
    Texture index = mp.getFloatTexture("index", 1.5);
    Texture bumpMap = mp.getFloatTextureOrNull("bumpmap");
    return new GlassMaterial(Kr, Kt, index, bumpMap);
  }

  BSDF getBSDF(DifferentialGeometry dgGeom, DifferentialGeometry dgShading) {
    DifferentialGeometry dgs;

    if (bumpMap != null) {
      dgs = new DifferentialGeometry();
      Material.Bump(bumpMap, dgGeom, dgShading, dgs);
    } else {
      dgs = dgShading;
    }

    double ior = index.evaluate(dgs);
    BSDF bsdf = new BSDF(dgs, dgGeom.nn, ior);
    Spectrum R = Kr.evaluate(dgs).clamp();
    Spectrum T = Kt.evaluate(dgs).clamp();

    if (!R.isBlack()) {
        bsdf.add(new SpecularReflection(R,
                    new FresnelDielectric(1.0, ior)));
    }

    if (!T.isBlack()) {
      bsdf.add(new SpecularTransmission(T, 1.0, ior));
    }

    return bsdf;
  }

  Texture Kr;
  Texture Kt;
  Texture index;
  Texture bumpMap;
}
