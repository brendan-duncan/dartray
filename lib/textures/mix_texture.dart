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
part of textures;

class MixTexture extends Texture {
  MixTexture(this.tex1, this.tex2, this.amount);

  static MixTexture CreateFloat(Transform tex2world, TextureParams tp) {
    return new MixTexture(
      tp.getFloatTexture("tex1", 0.0),
      tp.getFloatTexture("tex2", 1.0),
      tp.getFloatTexture("amount", 0.5));
  }

  static MixTexture CreateSpectrum(Transform tex2world, TextureParams tp) {
    return new MixTexture(
            tp.getSpectrumTexture("tex1", new Spectrum(0.0)),
            tp.getSpectrumTexture("tex2", new Spectrum(1.0)),
            tp.getFloatTexture("amount", 0.5));
  }

  evaluate(DifferentialGeometry dg) {
    var t1 = tex1.evaluate(dg);
    var t2 = tex2.evaluate(dg);
    double amt = amount.evaluate(dg);
    return t1 * (1.0 - amt) + t2 * amt;
  }

  Texture tex1, tex2;
  Texture amount;
}
