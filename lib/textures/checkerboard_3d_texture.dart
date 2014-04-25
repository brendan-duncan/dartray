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

class Checkerboard3DTexture extends Texture {
  Checkerboard3DTexture(this.mapping, this.tex1, this.tex2);

  evaluate(DifferentialGeometry dg) {
    Vector dpdx = new Vector();
    Vector dpdy = new Vector();
    Point p = mapping.map(dg, dpdx, dpdy);
    if ((p.x.floor() + p.y.floor() + p.z.floor()) % 2 == 0) {
      return tex1.evaluate(dg);
    } else {
      return tex2.evaluate(dg);
    }
  }

  TextureMapping3D mapping;
  Texture tex1, tex2;
}
