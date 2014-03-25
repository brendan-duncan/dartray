/****************************************************************************
 *  Copyright (C) 2014 by Brendan Duncan.                                   *
 *                                                                          *
 *  This file is part of DartRay.                                           *
 *                                                                          *
 *  Licensed under the Apache License, Version 2.0 (the "License");         *
 *  you may not use this file except in compliance with the License.        *
 *  You may obtain a copy of the License at                                 *
 *                                                                          *
 *  http://www.apache.org/licenses/LICENSE-2.0                              *
 *                                                                          *
 *  Unless required by applicable law or agreed to in writing, software     *
 *  distributed under the License is distributed on an "AS IS" BASIS,       *
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 *  See the License for the specific language governing permissions and     *
 *  limitations under the License.                                          *
 *                                                                          *
 *   This project is based on PBRT v2 ; see http://www.pbrt.org             *
 *   pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.*
 ****************************************************************************/
part of core;

abstract class Primitive {
  Primitive() :
    primitiveId = _nextprimitiveId++;

  BBox worldBound();

  bool canIntersect() {
    return true;
  }

  bool intersect(Ray r, Intersection intersect);

  bool intersectP(Ray r);

  void refine(List<Primitive> refined) {
    LogSevere("Unimplemented Primitive.refine() method called!");
  }

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

  AreaLight getAreaLight();

  BSDF getBSDF(DifferentialGeometry dg, Transform ObjectToWorld);

  BSSRDF getBSSRDF(DifferentialGeometry dg, Transform ObjectToWorld);

  final int primitiveId;

  static int _nextprimitiveId = 1;
}
