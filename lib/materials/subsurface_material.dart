/****************************************************************************
 * Copyright (C) 2014 by Brendan Duncan.                                    *
 *                                                                          *
 * This file is part of DartRay.                                            *
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the 'License');          *
 * you may not use this file except in compliance with the License.         *
 * You may obtain a copy of the License at                                  *
 *                                                                          *
 * http://www.apache.org/licenses/LICENSE-2.0                               *
 *                                                                          *
 * Unless required by applicable law or agreed to in writing, software      *
 * distributed under the License is distributed on an 'AS IS' BASIS,        *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 * See the License for the specific language governing permissions and      *
 * limitations under the License.                                           *
 *                                                                          *
 * This project is based on PBRT v2 ; see http://www.pbrt.org               *
 * pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.  *
 ****************************************************************************/
part of materials;

class SubsurfaceMaterial extends Material {
  SubsurfaceMaterial(this.scale, this.Kr, this.sigma_a, this.sigma_prime_s,
                     this.eta, this.bumpMap);

  static SubsurfaceMaterial Create(Transform xform, TextureParams mp) {
    List<double> sa_rgb = [0.0011, 0.0024, 0.014];
    List<double> sps_rgb = [2.55, 3.21, 3.77];
    Spectrum sa = new Spectrum.rgb(sa_rgb[0], sa_rgb[1], sa_rgb[2]);
    Spectrum sps = new Spectrum.rgb(sps_rgb[0], sps_rgb[1], sps_rgb[2]);
    String name = mp.findString('name');
    bool found = GetVolumeScatteringProperties(name, sa, sps);
    if (name != '' && !found) {
      LogWarning('Named material \'$name\' not found.  Using defaults.');
    }

    double scale = mp.findFloat('scale', 1.0);
    Texture sigma_a = mp.getSpectrumTexture('sigma_a', sa);
    Texture sigma_prime_s = mp.getSpectrumTexture('sigma_prime_s', sps);
    Texture ior = mp.getFloatTexture('index', 1.3);
    Texture Kr = mp.getSpectrumTexture('Kr', new Spectrum(1.0));
    Texture bumpMap = mp.getFloatTextureOrNull('bumpmap');
    return new SubsurfaceMaterial(scale, Kr, sigma_a, sigma_prime_s, ior,
        bumpMap);
  }

  BSDF getBSDF(DifferentialGeometry dgGeom, DifferentialGeometry dgShading) {
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
    return new BSSRDF(sigma_a.evaluate(dgShading) * scale,
        sigma_prime_s.evaluate(dgShading) * scale, e);
  }

  double scale;
  Texture Kr;
  Texture sigma_a;
  Texture sigma_prime_s;
  Texture eta;
  Texture bumpMap;
}
