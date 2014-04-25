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

class Sample extends CameraSample {
  List<int> n1D = [];
  List<int> n2D = [];
  List<Float32List> oneD = [];
  List<Float32List> twoD = [];

  Sample(Sampler sampler, SurfaceIntegrator surf, VolumeIntegrator vol,
         Scene scene) {
    if (surf != null) {
      surf.requestSamples(sampler, this, scene);
    }
    if (vol != null) {
      vol.requestSamples(sampler, this, scene);
    }
    _allocateSampleMemory();
  }

  Sample.from(Sample other) {
    n1D = new List<int>.from(other.n1D);
    n2D = new List<int>.from(other.n2D);
    _allocateSampleMemory();
  }

  int add1D(int num) {
    n1D.add(num);
    return n1D.length - 1;
  }

  int add2D(int num) {
    n2D.add(num);
    return n2D.length - 1;
  }

  List<Sample> duplicate(int count) {
    List<Sample> ret = new List<Sample>(count);
    for (int i = 0; i < count; ++i) {
      ret[i] = new Sample.from(this);
    }
    return ret;
  }

  void _allocateSampleMemory() {
    oneD = new List<Float32List>(n1D.length);
    twoD = new List<Float32List>(n2D.length);

    // Compute total number of sample values needed
    int totSamples = 0;

    for (int i = 0; i < n1D.length; ++i) {
      oneD[i] = new Float32List(n1D[i]);
    }

    for (int i = 0; i < n2D.length; ++i) {
      twoD[i] = new Float32List(2 * n2D[i]);
    }
  }
}
