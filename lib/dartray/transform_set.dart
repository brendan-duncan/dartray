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
part of dartray;

class TransformSet {
  TransformSet() {
    for (int i = 0; i < DartRay.MAX_TRANSFORMS; ++i) {
      t[i] = new Transform();
    }
  }

  TransformSet.from(TransformSet other) {
    t[0] = new Transform.from(other.t[0]);
    t[1] = new Transform.from(other.t[1]);
  }

  Transform operator [](int i) {
    assert(i >= 0 && i < DartRay.MAX_TRANSFORMS);
    return t[i];
  }

  operator []=(int i, Transform transform) {
    assert(i >= 0 && i < DartRay.MAX_TRANSFORMS);
    t[i] = transform;
  }

  bool isAnimated() {
    for (int i = 0; i < DartRay.MAX_TRANSFORMS - 1; ++i) {
      if (t[i] != t[i + 1]) {
        return true;
      }
    }
    return false;
  }

  static TransformSet Inverse(TransformSet ts) {
    TransformSet t2 = new TransformSet();
    for (int i = 0; i < DartRay.MAX_TRANSFORMS; ++i) {
      t2.t[i] = Transform.Inverse(ts.t[i]);
    }
    return t2;
  }

  List<Transform> t = new List<Transform>(DartRay.MAX_TRANSFORMS);
}
