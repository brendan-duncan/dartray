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

class ShinyMetalMaterial extends Material {
  ShinyMetalMaterial(this.Ks, this.roughness, this.Kr, this.bumpMap);

  static ShinyMetalMaterial Create(Transform xform, TextureParams mp) {
    Texture Kr = mp.getSpectrumTexture('Kr', new Spectrum(1.0));
    Texture Ks = mp.getSpectrumTexture('Ks', new Spectrum(1.0));
    Texture roughness = mp.getFloatTexture('roughness', 0.1);
    Texture bumpMap = mp.getFloatTextureOrNull('bumpmap');

    return new ShinyMetalMaterial(Ks, roughness, Kr, bumpMap);
  }

  BSDF getBSDF(DifferentialGeometry dgGeom,
               DifferentialGeometry dgShading) {
    // Allocate _BSDF_, possibly doing bump-mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap != null) {
      dgs = new DifferentialGeometry();
      Material.Bump(bumpMap, dgGeom, dgShading, dgs);
    } else {
      dgs = dgShading;
    }

    BSDF bsdf = new BSDF(dgs, dgGeom.nn);
    Spectrum spec = Ks.evaluate(dgs).clamp();
    double rough = roughness.evaluate(dgs);
    Spectrum R = Kr.evaluate(dgs).clamp();

    MicrofacetDistribution md = new Blinn(1.0 / rough);
    Spectrum k = new Spectrum(0.0);

    if (!spec.isBlack()) {
      Fresnel frMf = new FresnelConductor(FresnelApproxEta(spec), k);
      bsdf.add(new Microfacet(new Spectrum(1.0), frMf, md));
    }

    if (!R.isBlack()) {
      Fresnel frSr = new FresnelConductor(FresnelApproxEta(R), k);
      bsdf.add(new SpecularReflection(new Spectrum(1.0), frSr));
    }

    return bsdf;
  }

  static Spectrum FresnelApproxEta(Spectrum Fr) {
    Spectrum reflectance = Fr.clamp(0.0, 0.999);
    Spectrum sqrtRefl = reflectance.sqrt();
    Spectrum one = new Spectrum(1.0);
    return (one + sqrtRefl) / (one - sqrtRefl);
  }

  static Spectrum FresnelApproxK(Spectrum Fr) {
    Spectrum reflectance = Fr.clamp(0.0, 0.999);
    Spectrum r = reflectance / (new Spectrum(1.0) - reflectance);
    return r.sqrt() * 2.0;
  }

  Texture Ks;
  Texture Kr;
  Texture roughness;
  Texture bumpMap;
}
