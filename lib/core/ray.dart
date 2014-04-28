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
 * A ray is a line with one end, starting at an [origin] and going off in
 * a [direction] from that [origin].
 */
class Ray {
  /// The starting point of the ray.
  Point origin;
  /// The direction the ray is pointing toward.
  Vector direction;
  /// The minimum distance for a valid intersection test with this ray.
  /// This is used to reduce self-intersections with geometry.
  double minDistance;
  /// The maximum distance for a valid intersection test.
  double maxDistance;
  /// The animation time of this ray, for testing against animated scenes.
  double time;
  /// The recursive depth of this ray, for secondary rays.
  int depth;

  Ray([Point o, Vector d, this.minDistance = 0.0, this.maxDistance = INFINITY,
      this.time = 0.0, this.depth = 0])
      : origin = o == null ? new Point() : new Point.from(o),
        direction = d == null ? new Vector() : new Vector.from(d);

  Ray.from(Ray other)
      : origin = new Point.from(other.origin),
        direction = new Vector.from(other.direction),
        minDistance = other.minDistance,
        maxDistance = other.maxDistance,
        time = other.time,
        depth = other.depth;

  Ray.withParent(Point origin, Vector direction, Ray parent,
                 [this.minDistance = 0.0, this.maxDistance = INFINITY])
      : this.origin = new Point.from(origin),
        this.direction = new Vector.from(direction),
        this.time = parent.time,
        this.depth = parent.depth + 1;

  void set(Point o, Vector d, double minDist, double maxDist, double time) {
    origin = new Point.from(o);
    direction = new Vector.from(d);
    minDistance = minDist;
    maxDistance = maxDist;
    this.time = time;
  }

  Point pointAt(double t) =>
      new Point.from(origin + (direction * t));

  bool hasNaNs() => origin.hasNaNs() || direction.hasNaNs() ||
                    minDistance.isNaN || maxDistance.isNaN;
}
