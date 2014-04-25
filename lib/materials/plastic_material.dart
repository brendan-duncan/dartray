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

class PlasticMaterial extends Material {
  PlasticMaterial(this.Kd, this.Ks, this.roughness, this.bumpMap);

  static PlasticMaterial Create(Transform xform, TextureParams mp) {
    Texture Kd = mp.getSpectrumTexture("Kd", new Spectrum(0.25));
    Texture Ks = mp.getSpectrumTexture("Ks", new Spectrum(0.25));
    Texture roughness = mp.getFloatTexture("roughness", 0.1);
    Texture bumpMap = mp.getFloatTextureOrNull("bumpmap");
    return new PlasticMaterial(Kd, Ks, roughness, bumpMap);
  }

  BSDF getBSDF(DifferentialGeometry dgGeom,
               DifferentialGeometry dgShading) {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs = new DifferentialGeometry();
    if (bumpMap != null) {
      Material.Bump(bumpMap, dgGeom, dgShading, dgs);
    } else {
      dgs = dgShading;
    }

    BSDF bsdf = new BSDF(dgs, dgGeom.nn);
    Spectrum kd = Kd.evaluate(dgs).clamp();
    if (!kd.isBlack()) {
      BxDF diff = new Lambertian(kd);
      bsdf.add(diff);
    }

    Spectrum ks = Ks.evaluate(dgs).clamp();
    if (!ks.isBlack()) {
      Fresnel fresnel = new FresnelDielectric(1.5, 1.0);
      double rough = roughness.evaluate(dgs);
      BxDF spec = new Microfacet(ks, fresnel, new Blinn(1.0 / rough));
      bsdf.add(spec);
    }

    return bsdf;
  }

  Texture Kd;
  Texture Ks;
  Texture roughness;
  Texture bumpMap;
}
