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

class ShapeSet {
  ShapeSet(Shape s) {
    List<Shape> todo = [];
    todo.add(s);
    while (todo.isNotEmpty) {
      Shape sh = todo.last;
      todo.removeLast();
      if (sh.canIntersect()) {
        shapes.add(sh);
      } else {
        sh.refine(todo);
      }
    }

    if (shapes.length > 64) {
      LogWarning('Area light geometry turned into ${shapes.length} shapes; '
                 'may be very inefficient.');
    }

    // Compute total area of shapes in _ShapeSet_ and area CDF
    area = 0.0;
    for (int i = 0; i < shapes.length; ++i) {
      double a = shapes[i].area();
      areas.add(a);
      area += a;
    }

    areaDistribution = new Distribution1D(areas, areas.length);
  }

  Point sample(LightSample ls, Normal Ns, [Point p]) {
    if (p == null) {
      int sn = areaDistribution.sampleDiscrete(ls.uComponent, null) %
               shapes.length;
      return shapes[sn].sample(ls.uPos[0], ls.uPos[1], Ns);
    }

    int sn = areaDistribution.sampleDiscrete(ls.uComponent, null) %
             shapes.length;

    Point pt = shapes[sn].sample2(p, ls.uPos[0], ls.uPos[1], Ns);
    // Find closest intersection of ray with shapes in ShapeSet
    Ray r = new Ray(p, pt - p, 1.0e-3, INFINITY);
    List<double> rayEps = [0.0];
    List<double> thit = [1.0];
    bool anyHit = false;
    DifferentialGeometry dg = new DifferentialGeometry();

    for (int i = 0; i < shapes.length; ++i) {
      anyHit = shapes[i].intersect(r, thit, rayEps, dg) || anyHit;
    }

    if (anyHit) {
      Ns.copy(dg.nn);
    }

    return r.pointAt(thit[0]);
  }

  double pdf(Point p, [Vector wi]) {
    if (wi != null) {
      double pdf = 0.0;
      for (int i = 0; i < shapes.length; ++i) {
        pdf += areas[i] * shapes[i].pdf2(p, wi);
      }
      return pdf / area;
    }

    double pdf = 0.0;
    for (int i = 0; i < shapes.length; ++i) {
      pdf += areas[i] * shapes[i].pdf(p);
    }
    return pdf / area;
  }

  List<Shape> shapes = [];
  double area;
  List<double> areas = [];
  Distribution1D areaDistribution;
}
