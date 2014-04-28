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
 * The base class for geometric shapes.
 */
abstract class Shape {
  final Transform objectToWorld;
  final Transform worldToObject;
  bool reverseOrientation;
  bool transformSwapsHandedness = false;
  final int shapeId;

  Shape(this.objectToWorld, this.worldToObject, this.reverseOrientation)
      : shapeId = _nextShapeId++;

  BBox objectBound();

  BBox worldBound() {
    return objectToWorld.transformBBox(objectBound());
  }

  bool canIntersect() {
    return true;
  }

  void refine(List<Shape> refined) {
    LogSevere('Unimplemented Shape.refine() method called');
  }

  /**
   * Tests for an intersection of the given [ray] with this geometry,
   * returning true if an intersection was found.  If the intersection was
   * found, [tHit], [rayEpsilon] and [dg] will be set to information about
   * the intersection point.
   */
  bool intersect(Ray ray, List<double> tHit, List<double> rayEpsilon,
                 DifferentialGeometry dg) {
    LogSevere('Unimplemented Shape.intersect() method called');
    return false;
  }

  /**
   * Tests for an intersection of the given [ray] with this geometry, without
   * computing any information about the actual intersection.  Because it
   * doesn't compute the intersection data, it is more efficient than
   * [intersect].
   */
  bool intersectP(Ray ray) {
    LogSevere('Unimplemented Shape.intersectP() method called');
    return false;
  }

  void getShadingGeometry(Transform obj2world,
                          DifferentialGeometry dg,
                          DifferentialGeometry dgShading) {
    dgShading.copy(dg);
  }

  /**
   * Surface area of the geometry.
   */
  double area() {
    LogSevere('Unimplemented Shape.area() method called');
    return 0.0;
  }

  Point sample(double u1, double u2, Normal Ns) {
    LogSevere('Unimplemented Shape::sample() method called');
    return new Point();
  }

  double pdf(Point Pshape) {
    return 1.0 / area();
  }

  Point sample2(Point P, double u1, double u2, Normal Ns) {
    return sample(u1, u2, Ns);
  }

  double pdf2(Point p, Vector wi) {
    // Intersect sample ray with area light geometry
    DifferentialGeometry dgLight = new DifferentialGeometry();
    Ray ray = new Ray(p, wi, 1.0e-3);
    ray.depth = -1; // temporary hack to ignore alpha mask

    List<double> thit = [0.0];
    List<double> rayEpsilon = [0.0];
    if (!intersect(ray, thit, rayEpsilon, dgLight)) {
      return 0.0;
    }

    // Convert light sample weight to solid angle measure
    double pdf = Vector.DistanceSquared(p, ray.pointAt(thit[0])) /
                 (Vector.AbsDot(dgLight.nn, -wi) * area());

    if (pdf.isInfinite) {
      pdf = 0.0;
    }

    return pdf;
  }

  static int _nextShapeId = 1;
}
