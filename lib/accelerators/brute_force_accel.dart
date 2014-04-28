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
part of accelerators;

class BruteForceAccel extends Aggregate {
  BruteForceAccel(List<Primitive> p) {
    LogInfo('Building Brute Force Acceleration Structures.');
    for (int i = 0; i < p.length; ++i) {
      p[i].fullyRefine(primitives);
    }

    // Compute bounds and choose grid resolution
    for (int i = 0; i < primitives.length; ++i) {
      bounds = BBox.Union(bounds, primitives[i].worldBound());
    }
  }

  BBox worldBound() => bounds;

  bool canIntersect() {
    return true;
  }

  bool intersect(Ray ray, Intersection isect) {
    List<double> rayT = [0.0];
    if (bounds.inside(ray.pointAt(ray.minDistance))) {
      rayT[0] = ray.minDistance;
    } else if (!bounds.intersectP(ray, rayT)) {
      return false;
    }

    bool hitSomething = false;
    for (int i = 0; i < primitives.length; ++i) {
      if (primitives[i].intersect(ray, isect)) {
        hitSomething = true;
      }
    }
    return hitSomething;
  }

  bool intersectP(Ray ray) {
    for (int i = 0; i < primitives.length; ++i) {
      if (primitives[i].intersectP(ray)) {
        return true;
      }
    }

    return false;
  }

  static BruteForceAccel Create(List<Primitive> prims, ParamSet ps) {
    return new BruteForceAccel(prims);
  }

  List<Primitive> primitives = [];
  BBox bounds = new BBox();
}
