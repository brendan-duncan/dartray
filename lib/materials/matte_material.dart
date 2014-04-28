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

class MatteMaterial extends Material {
  MatteMaterial(this.Kd, this.sigma, this.bumpMap);

  static MatteMaterial Create(Transform xform, TextureParams mp) {
    Texture Kd = mp.getSpectrumTexture('Kd', new Spectrum(0.5));
    Texture sigma = mp.getFloatTexture('sigma', 0.0);
    Texture bumpMap = mp.getFloatTextureOrNull('bumpmap');
    return new MatteMaterial(Kd, sigma, bumpMap);
  }

  BSDF getBSDF(DifferentialGeometry dgGeom,
               DifferentialGeometry dgShading) {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap != null) {
      dgs = new DifferentialGeometry();
      Material.Bump(bumpMap, dgGeom, dgShading, dgs);
    } else {
      dgs = dgShading;
    }

    BSDF bsdf = new BSDF(dgs, dgGeom.nn);

    // Evaluate textures for _MatteMaterial_ material and allocate BRDF
    Spectrum r = Kd.evaluate(dgs).clamp();
    double sig = sigma.evaluate(dgs).clamp(0.0, 90.0);
    if (!r.isBlack()) {
      if (sig == 0.0) {
        bsdf.add(new Lambertian(r));
      } else {
        bsdf.add(new OrenNayar(r, sig));
      }
    }

    return bsdf;
  }

  Texture Kd; // Texture<RGBSpectrum>
  Texture sigma; // Texture<double>
  Texture bumpMap; // Texture<double>
}
