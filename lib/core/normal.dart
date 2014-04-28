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

/**
 * A vector that's interpreted and transformed as a surface normal.
 */
class Normal extends Vector {
  static final Normal ZERO = new Normal(0.0, 0.0, 0.0);

  Normal([num x = 0.0, num y = 0.0, num z = 0.0])
      : super(x, y, z);

  Normal.from(Vector other)
      : super.from(other);

  Normal normalize() {
    invScale(length());
    return this;
  }

  Normal operator*(num s) =>
      new Normal(x * s, y * s, z * s);

  Normal operator/(num s) =>
      new Normal(x / s, y / s, z / s);

  Normal operator+(Vector p) =>
      new Normal(x + p.x, y + p.y, z + p.z);

  Normal operator-(Vector p) =>
      new Normal(x - p.x, y - p.y, z - p.z);

  Normal operator-() =>
        new Normal(-x, -y, -z);

  static Normal Normalize(Vector n) {
    return new Normal.from(n).normalize();
  }

  static Normal FaceForward(Normal n, Vector n2) {
    return (Vector.Dot(n, n2) < 0.0) ? new Normal(-n.x, -n.y, -n.z) : n;
  }
}

