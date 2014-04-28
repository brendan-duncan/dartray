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
 * Stores information about the point on a surface from a ray intersection.
 * This includes the point, normal, and various derivatives of the surface.
 */
class DifferentialGeometry {
  Point p;
  Normal nn;
  double u;
  double v;
  Shape shape;
  Vector dpdu;
  Vector dpdv;
  Normal dndu;
  Normal dndv;
  Vector dpdx;
  Vector dpdy;
  double dudx;
  double dvdx;
  double dudy;
  double dvdy;

  DifferentialGeometry()
      : p = new Point(),
        nn = new Normal(),
        u = 0.0,
        v = 0.0,
        dpdu = new Vector(),
        dpdv = new Vector(),
        dndu = new Normal(),
        dndv = new Normal(),
        dpdx = new Vector(),
        dpdy = new Vector(),
        dudx = 0.0,
        dvdx = 0.0,
        dudy = 0.0,
        dvdy = 0.0;

  DifferentialGeometry.from(DifferentialGeometry other)
      : p = new Point.from(other.p),
        nn = new Normal.from(other.nn),
        u = other.u,
        v = other.v,
        shape = other.shape,
        dpdu = new Vector.from(other.dpdu),
        dpdv = new Vector.from(other.dpdv),
        dndu = new Normal.from(other.dndu),
        dndv = new Normal.from(other.dndv),
        dpdx = new Vector.from(other.dpdx),
        dpdy = new Vector.from(other.dpdy),
        dudx = other.dudx,
        dvdx = other.dvdx,
        dudy = other.dudy,
        dvdy = other.dvdy;

  DifferentialGeometry set(Point p, Vector dpdu, Vector dpdv, Normal dndu,
                           Normal dndv, double u, double v, Shape shape) {
    this.p = p;
    this.dpdu = dpdu;
    this.dpdv = dpdv;
    this.dndu = dndu;
    this.dndv = dndv;
    this.nn = Vector.Normalize(Vector.Cross(dpdu, dpdv));
    this.u = u;
    this.v = v;
    this.shape = shape;
    dudx = 0.0;
    dvdx = 0.0;
    dudy = 0.0;
    dvdy = 0.0;

    // Adjust normal based on orientation and handedness
    int a = shape != null && shape.reverseOrientation ? 1 : 0;
    int b = shape != null && shape.transformSwapsHandedness ? 1 : 0;

    if (shape != null && (a ^ b) != 0) {
      nn *= -1.0;
    }

    return this;
  }

  void copy(DifferentialGeometry other) {
    p.copy(other.p);
    nn.copy(other.nn);
    u = other.u;
    v = other.v;
    shape = other.shape;
    dpdu.copy(other.dpdu);
    dpdv.copy(other.dpdv);
    dndu.copy(other.dndu);
    dndv.copy(other.dndv);
    dpdx.copy(other.dpdx);
    dpdy.copy(other.dpdy);
    dudx = other.dudx;
    dvdx = other.dvdx;
    dudy = other.dudy;
    dvdy = other.dvdy;
  }

  void computeDifferentials(RayDifferential ray) {
    if (ray.hasDifferentials) {
      // Estimate screen space change in pt and (u,v)

      // Compute auxiliary intersection points with plane
      double d = -Vector.Dot(nn, new Vector(p.x, p.y, p.z));
      Vector rxv = new Vector(ray.rxOrigin.x, ray.rxOrigin.y, ray.rxOrigin.z);
      double tx = -(Vector.Dot(nn, rxv) + d) / Vector.Dot(nn, ray.rxDirection);
      if (tx.isNaN) {
        dudx = dvdx = 0.0;
        dudy = dvdy = 0.0;
        dpdx = new Vector(0.0, 0.0, 0.0);
        dpdy = new Vector(0.0, 0.0, 0.0);
        return;
      }

      Point px = new Point.from(ray.rxOrigin + ray.rxDirection * tx);
      Vector ryv = new Vector.from(ray.ryOrigin);
      double ty = -(Vector.Dot(nn, ryv) + d) / Vector.Dot(nn, ray.ryDirection);
      if (ty.isNaN) {
        dudx = dvdx = 0.0;
        dudy = dvdy = 0.0;
        dpdx = new Vector(0.0, 0.0, 0.0);
        dpdy = new Vector(0.0, 0.0, 0.0);
        return;
      }

      Point py = new Point.from(ray.ryOrigin + ray.ryDirection * ty);
      dpdx = px - p;
      dpdy = py - p;

      // Compute (u,v) offsets at auxiliary points

      // Initialize A, Bx, and By matrices for offset computation
      List<double> A = new List<double>(4);
      List<double> Bx = new List<double>(2);
      List<double> By = new List<double>(2);

      List<int> axes = new List<int>(2);
      if (nn.x.abs() > nn.y.abs() && nn.x.abs() > nn.z.abs()) {
        axes[0] = 1;
        axes[1] = 2;
      } else if (nn.y.abs() > nn.z.abs()) {
        axes[0] = 0;
        axes[1] = 2;
      } else {
        axes[0] = 0;
        axes[1] = 1;
      }

      // Initialize matrices for chosen projection plane
      A[0] = dpdu[axes[0]];
      A[1] = dpdv[axes[0]];
      A[2] = dpdu[axes[1]];
      A[3] = dpdv[axes[1]];
      Bx[0] = px[axes[0]] - p[axes[0]];
      Bx[1] = px[axes[1]] - p[axes[1]];
      By[0] = py[axes[0]] - p[axes[0]];
      By[1] = py[axes[1]] - p[axes[1]];

      List<double> _du = [0.0];
      List<double> _dv = [0.0];
      if (!SolveLinearSystem2x2(A, Bx, _du, _dv)) {
        dudx = 0.0;
        dvdx = 0.0;
      } else {
        dudx = _du[0];
        dvdx = _dv[0];
      }

      if (!SolveLinearSystem2x2(A, By, _du, _dv)) {
        dudy = 0.0;
        dvdy = 0.0;
      } else {
        dudy = _du[0];
        dvdy = _dv[0];
      }
    } else {
      dudx = dvdx = 0.0;
      dudy = dvdy = 0.0;
      dpdx = new Vector(0.0, 0.0, 0.0);
      dpdy = new Vector(0.0, 0.0, 0.0);
    }
  }
}
