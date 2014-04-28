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

class BSDFSample {
  final Float32List uDir;
  double uComponent = 0.0;

  BSDFSample([double up0 = 0.0, double up1 = 0.0, this.uComponent = 0.0])
      : uDir = new Float32List(2) {
    uDir[0] = up0;
    uDir[1] = up1;
  }

  BSDFSample.from(BSDFSample other)
      : uDir = new Float32List.fromList(other.uDir),
        uComponent = other.uComponent;

  BSDFSample.random(RNG rng)
      : uDir = new Float32List(2) {
    uDir[0] = rng.randomFloat();
    uDir[1] = rng.randomFloat();
    uComponent = rng.randomFloat();
  }

  BSDFSample.sample(Sample sample, BSDFSampleOffsets offsets, int n)
      : uDir = new Float32List(2) {
    uDir[0] = sample.twoD[offsets.dirOffset][2 * n];
    uDir[1] = sample.twoD[offsets.dirOffset][2 * n + 1];
    uComponent = sample.oneD[offsets.componentOffset][n];
  }
}
