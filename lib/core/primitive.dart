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
 * Base class for objects that can be for ray intersections.
 */
abstract class Primitive {
  Primitive()
      : primitiveId = _nextprimitiveId++;

  /**
   * The world-space bounding box for this primitive.
   */
  BBox worldBound();

  /**
   * If false, the primitive should be refined using the [refine] method into
   * its sub primitives, until primitives that can be ray intersected are
   * found.
   */
  bool canIntersect() {
    return true;
  }

  /**
   * Test the primitive for an intersection with the given [ray]. This will
   * also compute the geometric data at the intersection point, such as
   * the normal and derivatives.
   */
  bool intersect(Ray ray, Intersection intersect);

  /**
   * A faster version of [intersect] that only tests for [ray] intersection,
   * but does not compute intersection geometric data.
   */
  bool intersectP(Ray ray);

  /**
   * If the primitive is not directly intersectable, then it should be
   * refined into sub-primitives. This method will refine the primitive one
   * level, which may contain non-intersectable sub-primitives. Use
   * [fullyRefine] to get the list of all intersectable sub-primitives.
   */
  void refine(List<Primitive> refined) {
    LogSevere("Unimplemented Primitive.refine() method called!");
  }

  /**
   * Recursively refine the primitive until all intersectable sub-primitives
   * are found.
   */
  void fullyRefine(List<Primitive> refined) {
    List<Primitive> todo = [];
    todo.add(this);
    while (todo.isNotEmpty) {
      // Refine last primitive in todo list
      Primitive prim = todo.last;
      todo.removeLast();
      if (prim.canIntersect()) {
        refined.add(prim);
      } else {
        prim.refine(todo);
      }
    }
  }

  /**
   * If the primitive is an area light, it will have an [AreaLight] associated
   * with it.
   */
  AreaLight getAreaLight();

  /**
   * Return the bi-directional scattering distribution function (BSDF) to
   * shade the surface of the primitive at the given point.
   */
  BSDF getBSDF(DifferentialGeometry dg, Transform ObjectToWorld);

  /**
   * Return the Bidirectional surface scattering reflectance distribution
   * function (BSSRDF) to compute things like sub-surface scattering.
   */
  BSSRDF getBSSRDF(DifferentialGeometry dg, Transform ObjectToWorld);

  /// A unique identifier for the primitive, useful for debugging.
  final int primitiveId;

  static int _nextprimitiveId = 1;
}
