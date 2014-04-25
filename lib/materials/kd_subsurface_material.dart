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

class KdSubsurfaceMaterial extends Material {
  KdSubsurfaceMaterial(this.Kd, this.Kr, this.meanfreepath, this.eta,
                       this.bumpMap);

  static KdSubsurfaceMaterial Create(Transform xform, TextureParams mp) {
    Texture kd = mp.getSpectrumTexture("Kd", new Spectrum(0.5));
    Texture mfp = mp.getFloatTexture("meanfreepath", 1.0);
    Texture ior = mp.getFloatTexture("index", 1.3);
    Texture kr = mp.getSpectrumTexture("Kr", new Spectrum(1.0));
    Texture bumpMap = mp.getFloatTextureOrNull("bumpmap");
    return new KdSubsurfaceMaterial(kd, kr, mfp, ior, bumpMap);
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
    Spectrum R = Kr.evaluate(dgs).clamp();
    double e = eta.evaluate(dgs);

    if (!R.isBlack()) {
      bsdf.add(new SpecularReflection(R, new FresnelDielectric(1.0, e)));
    }

    return bsdf;
  }

  BSSRDF getBSSRDF(DifferentialGeometry dgGeom,
                   DifferentialGeometry dgShading) {
    double e = eta.evaluate(dgShading);
    double mfp = meanfreepath.evaluate(dgShading);
    Spectrum kd = Kd.evaluate(dgShading).clamp();
    Spectrum sigma_a = new Spectrum(0.0);
    Spectrum sigma_prime_s = new Spectrum(0.0);
    SubsurfaceFromDiffuse(kd, mfp, e, sigma_a, sigma_prime_s);

    return new BSSRDF(sigma_a, sigma_prime_s, e);
  }

  Texture Kd, Kr;
  Texture meanfreepath, eta, bumpMap;
}
