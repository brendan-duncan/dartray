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

class TranslucentMaterial extends Material {
  TranslucentMaterial(this.Kd, this.Ks, this.roughness, this.reflect,
                      this.transmit, this.bumpMap);

  static TranslucentMaterial Create(Transform xform, TextureParams mp) {
    Texture Kd = mp.getSpectrumTexture("Kd", new Spectrum(0.25));
    Texture Ks = mp.getSpectrumTexture("Ks", new Spectrum(0.25));
    Texture reflect = mp.getSpectrumTexture("reflect", new Spectrum(0.5));
    Texture transmit = mp.getSpectrumTexture("transmit", new Spectrum(0.5));
    Texture roughness = mp.getFloatTexture("roughness", 0.1);
    Texture bumpMap = mp.getFloatTextureOrNull("bumpmap");
    return new TranslucentMaterial(Kd, Ks, roughness, reflect, transmit,
                                   bumpMap);
  }

  BSDF getBSDF(DifferentialGeometry dgGeom, DifferentialGeometry dgShading) {
    double ior = 1.5;
    DifferentialGeometry dgs = new DifferentialGeometry();

    if (bumpMap != null) {
      Material.Bump(bumpMap, dgGeom, dgShading, dgs);
    } else {
      dgs = dgShading;
    }

    BSDF bsdf = new BSDF(dgs, dgGeom.nn, ior);

    Spectrum r = reflect.evaluate(dgs).clamp();
    Spectrum t = transmit.evaluate(dgs).clamp();

    if (r.isBlack() && t.isBlack()) {
      return bsdf;
    }

    Spectrum kd = Kd.evaluate(dgs).clamp();

    if (!kd.isBlack()) {
      if (!r.isBlack()) {
        bsdf.add(new Lambertian(r * kd));
      }
      if (!t.isBlack()) {
        bsdf.add(new BRDFToBTDF(new Lambertian(t * kd)));
      }
    }

    Spectrum ks = Ks.evaluate(dgs).clamp();

    if (!ks.isBlack()) {
      double rough = roughness.evaluate(dgs);
      if (!r.isBlack()) {
        Fresnel fresnel = new FresnelDielectric(ior, 1.0);
        bsdf.add(new Microfacet(r * ks, fresnel, new Blinn(1.0 / rough)));
      }

      if (!t.isBlack()) {
        Fresnel fresnel = new FresnelDielectric(ior, 1.0);
        bsdf.add(new BRDFToBTDF(new Microfacet(t * ks, fresnel,
            new Blinn(1.0 / rough))));
      }
    }

    return bsdf;
  }

  Texture Kd;
  Texture Ks;
  Texture roughness;
  Texture reflect;
  Texture transmit;
  Texture bumpMap;
}
