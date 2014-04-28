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
part of core;

class LightSample {
  LightSample([double up0 = 0.0, double up1 = 0.0, this.uComponent = 0.0])
      : uPos = new Float32List(2) {
    uPos[0] = up0;
    uPos[1] = up1;
  }

  LightSample.from(LightSample other)
      : uPos = new Float32List.fromList(other.uPos),
        uComponent = other.uComponent;

  LightSample.sample(Sample sample, LightSampleOffsets offsets, int n)
      : uPos = new Float32List(2) {
    assert(n < sample.n2D[offsets.posOffset]);
    assert(n < sample.n1D[offsets.componentOffset]);
    uPos[0] = sample.twoD[offsets.posOffset][2 * n];
    uPos[1] = sample.twoD[offsets.posOffset][2 * n + 1];
    uComponent = sample.oneD[offsets.componentOffset][n];
    assert(uPos[0] >= 0.0 && uPos[0] <= 1.0);
    assert(uPos[1] >= 0.0 && uPos[1] <= 1.0);
    assert(uComponent >= 0.0 && uComponent <= 1.0);
  }

  LightSample.random(RNG rng)
      : uPos = new Float32List(2) {
    uPos[0] = rng.randomFloat();
    uPos[1] = rng.randomFloat();
    uComponent = rng.randomFloat();
  }

  final Float32List uPos;
  double uComponent;
}
