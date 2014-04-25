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
 * A primitive that has [Shape] geometry and a [Material]. All geometetry
 * primitives can also be made into an [AreaLight].
 */
class GeometricPrimitive extends Primitive {
  GeometricPrimitive(this.shape, this.material, this.areaLight);

  bool canIntersect() {
    return shape.canIntersect();
  }

  void refine(List<Primitive> refined) {
    List<Shape> r = [];
    shape.refine(r);
    for (int i = 0; i < r.length; ++i) {
      GeometricPrimitive gp = new GeometricPrimitive(r[i], material, areaLight);
      refined.add(gp);
    }
  }

  BBox worldBound() {
    return shape.worldBound();
  }

  bool intersect(Ray r, Intersection isect) {
    List<double> thit = [0.0];
    List<double> rayEpsilon = [0.0];
    if (!shape.intersect(r, thit, rayEpsilon, isect.dg)) {
      return false;
    }
    isect.primitive = this;
    isect.worldToObject = new Transform.from(shape.worldToObject);
    isect.objectToWorld = new Transform.from(shape.objectToWorld);
    isect.shapeId = shape.shapeId;
    isect.primitiveId = primitiveId;
    isect.rayEpsilon = rayEpsilon[0];
    r.maxDistance = thit[0];
    return true;
  }

  bool intersectP(Ray r) {
    return shape.intersectP(r);
  }

  AreaLight getAreaLight() {
    return areaLight;
  }

  BSDF getBSDF(DifferentialGeometry dg, Transform objectToWorld) {
    DifferentialGeometry dgs = new DifferentialGeometry();
    shape.getShadingGeometry(objectToWorld, dg, dgs);
    return material.getBSDF(dg, dgs);
  }

  BSSRDF getBSSRDF(DifferentialGeometry dg, Transform objectToWorld) {
    DifferentialGeometry dgs = new DifferentialGeometry();
    shape.getShadingGeometry(objectToWorld, dg, dgs);
    return material.getBSSRDF(dg, dgs);
  }

  Shape shape;
  Material material;
  AreaLight areaLight;
}
