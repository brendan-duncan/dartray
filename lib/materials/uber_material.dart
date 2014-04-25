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

class UberMaterial extends Material {
  UberMaterial(this.Kd, this.Ks, this.Kr, this.Kt, this.roughness,
               this.opacity, this.eta, this.bumpMap);

  static UberMaterial Create(Transform xform, TextureParams mp) {
    Texture Kd = mp.getSpectrumTexture("Kd", new Spectrum(0.25));
    Texture Ks = mp.getSpectrumTexture("Ks", new Spectrum(0.25));
    Texture Kr = mp.getSpectrumTexture("Kr", new Spectrum(0.0));
    Texture Kt = mp.getSpectrumTexture("Kt", new Spectrum(0.0));
    Texture roughness = mp.getFloatTexture("roughness", 0.1);
    Texture eta = mp.getFloatTexture("index", 1.5);
    Texture opacity = mp.getSpectrumTexture("opacity", new Spectrum(1.0));
    Texture bumpMap = mp.getFloatTextureOrNull("bumpmap");
    return new UberMaterial(Kd, Ks, Kr, Kt, roughness, opacity, eta, bumpMap);
  }

  BSDF getBSDF(DifferentialGeometry dgGeom, DifferentialGeometry dgShading) {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap != null) {
      dgs = new DifferentialGeometry();
      Material.Bump(bumpMap, dgGeom, dgShading, dgs);
    } else {
      dgs = dgShading;
    }

    BSDF bsdf = new BSDF(dgs, dgGeom.nn);

    Spectrum op = opacity.evaluate(dgs).clamp();
    if (!op.isValue(1.0)) {
      BxDF tr = new SpecularTransmission(-op + new Spectrum(1.0), 1.0, 1.0);
      bsdf.add(tr);
    }

    Spectrum kd = op * Kd.evaluate(dgs).clamp();
    if (!kd.isBlack()) {
      BxDF diff = new Lambertian(kd);
      bsdf.add(diff);
    }

    double e = eta.evaluate(dgs);
    Spectrum ks = op * Ks.evaluate(dgs).clamp();
    if (!ks.isBlack()) {
      Fresnel fresnel = new FresnelDielectric(e, 1.0);
      double rough = roughness.evaluate(dgs);
      BxDF spec = new Microfacet(ks, fresnel, new Blinn(1.0 / rough));
      bsdf.add(spec);
    }

    Spectrum kr = op * Kr.evaluate(dgs).clamp();
    if (!kr.isBlack()) {
      Fresnel fresnel = new FresnelDielectric(e, 1.0);
      bsdf.add(new SpecularReflection(kr, fresnel));
    }

    Spectrum kt = op * Kt.evaluate(dgs).clamp();
    if (!kt.isBlack()) {
      bsdf.add(new SpecularTransmission(kt, e, 1.0));
    }

    return bsdf;
  }

  Texture Kd, Ks, Kr, Kt, opacity;
  Texture roughness, eta, bumpMap;
}
