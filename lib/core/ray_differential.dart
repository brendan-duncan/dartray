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
 * A RayDifferential is a Ray that contains additional information used for
 * filtering textures.
 */
class RayDifferential extends Ray {
  bool hasDifferentials = false;
  Point rxOrigin = new Point();
  Point ryOrigin = new Point();
  Vector rxDirection = new Vector();
  Vector ryDirection = new Vector();

  RayDifferential([Point origin, Vector direction, double start = 0.0,
                   double end = INFINITY, double t = 0.0, int d = 0])
      : super(origin, direction, start, end, t, d);

  RayDifferential.fromRay(Ray ray)
      : super(new Point.from(ray.origin), new Vector.from(ray.direction),
              ray.minDistance, ray.maxDistance, ray.time, ray.depth);

  RayDifferential.child([Point origin, Vector direction, Ray parent,
                       double start = 0.0, end = INFINITY])
      : super.withParent(origin, direction, parent, start, end);

  void copy(RayDifferential b) {
    origin.copy(b.origin);
    direction.copy(b.direction);
    hasDifferentials = b.hasDifferentials;
    rxOrigin.copy(b.rxOrigin);
    ryOrigin.copy(b.ryOrigin);
    rxDirection.copy(b.rxDirection);
    ryDirection.copy(b.ryDirection);
  }

  void scaleDifferentials(double s) {
    rxOrigin = new Point.from(origin + (rxOrigin - origin) * s);
    ryOrigin = new Point.from(origin + (ryOrigin - origin) * s);
    rxDirection = direction + (rxDirection - direction) * s;
    ryDirection = direction + (ryDirection - direction) * s;
  }

  bool hasNaNs() {
    return super.hasNaNs() ||
         (hasDifferentials && (rxOrigin.hasNaNs() || ryOrigin.hasNaNs() ||
                               rxDirection.hasNaNs() || ryDirection.hasNaNs()));
  }
}
