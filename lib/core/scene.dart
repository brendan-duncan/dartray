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
 * A scene encompasses all of the geometry and lights to be rendered.
 */
class Scene {
  /// An aggregate primitive is one that contains sub-primitives, usually an
  /// acceleration structure.
  Primitive aggregate;
  /// The list of lights in the scene.
  List<Light> lights;
  /// A region for volumetric effects.
  VolumeRegion volumeRegion;
  /// The world-space boudning box of all geometry.
  BBox worldBound;

  Scene(this.aggregate, List<Light> lights, this.volumeRegion)
      : lights = new List<Light>.from(lights) {
    worldBound = aggregate.worldBound();
    if (volumeRegion != null) {
      worldBound = BBox.Union(worldBound, volumeRegion.worldBound());
    }
  }

  /**
   * Traces the given [ray] into the scene and returns true if the ray
   * intersected any of the primitives. If an intersection was found,
   * [isect] is filled in with information about the closest intersection
   * point along the ray.
   */
  bool intersect(Ray ray, Intersection isect) {
    Stats.STARTED_RAY_INTERSECTION(ray);
    bool hit = aggregate.intersect(ray, isect);
    Stats.FINISHED_RAY_INTERSECTION(ray, isect, hit);
    return hit;
  }

  /**
   * Traces the given [ray] into the scene and returns true if an intersection
   * was found.  Unlike [intersect], this does not return any information about
   * the intersection, and so is more efficient than [intersect].
   */
  bool intersectP(Ray ray) {
    Stats.STARTED_RAY_INTERSECTIONP(ray);
    bool hit = aggregate.intersectP(ray);
    Stats.FINISHED_RAY_INTERSECTIONP(ray, hit);
    return hit;
  }
}
