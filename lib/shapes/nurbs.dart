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
part of shapes;

class Nurbs extends Shape {
  Nurbs(Transform o2w, Transform w2o, bool reverseOrientation,
             this.nu, this.uorder, this.uknot, this.umin, this.umax,
             this.nv, this.vorder, this.vknot, this.vmin, this.vmax,
             this.P, this.isHomogeneous) :
    super(o2w, w2o, reverseOrientation);

  BBox objectBound() {
    if (!isHomogeneous) {
      // Compute object-space bound of non-homogeneous NURBS
      BBox bound = new BBox();
      for (int i = 0, pi = 0; i < nu * nv; ++i, pi += 3) {
        bound = BBox.UnionPoint(bound, new Point(P[pi], P[pi + 1], P[pi + 2]));
      }
      return bound;
    } else {
      // Compute object-space bound of homogeneous NURBS
      BBox bound = new BBox();
      for (int i = 0, pi = 0; i < nu * nv; ++i, pi += 4) {
        bound = BBox.UnionPoint(bound, new Point(P[pi] / P[pi + 3],
                                            P[pi + 1] / P[pi + 3],
                                            P[pi + 2] / P[pi + 3]));
      }
      return bound;
    }
  }

  BBox worldBound() {
    if (!isHomogeneous) {
      // Compute object-space bound of non-homogeneous NURBS
      BBox bound = new BBox();
      for (int i = 0, pi = 0; i < nu * nv; ++i, pi += 3) {
        Point pt = new Point(P[pi], P[pi + 1], P[pi + 2]);
        bound = BBox.UnionPoint(bound, objectToWorld.transformPoint(pt));
      }
      return bound;
    } else {
      // Compute object-space bound of homogeneous NURBS
      BBox bound = new BBox();
      for (int i = 0, pi = 0; i < nu * nv; ++i, pi += 4) {
        Point pt = new Point(P[pi] / P[pi + 3],
                             P[pi + 1] / P[pi + 3],
                             P[pi + 2] / P[pi + 3]);
        bound = BBox.UnionPoint(bound, objectToWorld.transformPoint(pt));
      }
      return bound;
    }
  }

  bool canIntersect() {
    return false;
  }

  void refine(List<Shape> refined) {
    // Compute NURBS dicing rates
    const int diceu = 30;
    const int dicev = 30;
    Float32List ueval = new Float32List(diceu);
    Float32List veval = new Float32List(dicev);
    List<Point> evalPs = new List<Point>(diceu * dicev);
    List<Normal> evalNs = new List<Normal>(diceu * dicev);

    for (int i = 0; i < diceu; ++i) {
      ueval[i] = Lerp(i / (diceu - 1), umin, umax);
    }

    for (int i = 0; i < dicev; ++i) {
      veval[i] = Lerp(i / (dicev - 1), vmin, vmax);
    }

    // Evaluate NURBS over grid of points
    Float32List uvs = new Float32List(2 * diceu * dicev);

    // Turn NURBS into triangles
    List<double> Pw = P;

    if (!isHomogeneous) {
      Pw = new Float32List(nu * nv * 4);
      for (int i = 0, pi = 0, wi = 0; i < nu * nv; ++i) {
        Pw[wi++] = P[pi++];
        Pw[wi++] = P[pi++];
        Pw[wi++] = P[pi++];
        Pw[wi++] = 1.0;
      }
    }

    for (int v = 0, pi = 0; v < dicev; ++v) {
      for (int u = 0; u < diceu; ++u, pi++) {
        uvs[2 * pi] = ueval[u];
        uvs[2 * pi + 1] = veval[v];

        Vector dPdu = new Vector();
        Vector dPdv = new Vector();
        Point pt = NurbsEvaluateSurface(uorder, uknot, nu, ueval[u],
                                        vorder, vknot, nv, veval[v],
                                        Pw, dPdu, dPdv);
        evalPs[pi] = pt;
        evalNs[pi] = Normal.Normalize(Vector.Cross(dPdu, dPdv));
      }
    }

    // Generate points-polygons mesh
    int nTris = 2 * (diceu - 1) * (dicev - 1);
    List<int> vertices = new Int32List(3 * nTris);
    int vi = 0;

    int VN(int u, int v) => v * diceu + u;

    // Compute the vertex offset numbers for the triangles
    for (int v = 0; v < dicev - 1; ++v) {
      for (int u = 0; u < diceu - 1; ++u) {
        vertices[vi++] = VN(u, v);
        vertices[vi++] = VN(u + 1, v);
        vertices[vi++] = VN(u + 1, v + 1);

        vertices[vi++] = VN(u, v);
        vertices[vi++] = VN(u + 1, v + 1);
        vertices[vi++] = VN(u, v + 1);
      }
    }

    int nVerts = diceu * dicev;

    ParamSet paramSet = new ParamSet();
    paramSet.addInt('indices', vertices);
    paramSet.addPoint('P', evalPs);
    paramSet.addFloat('uv', uvs);
    paramSet.addNormal('N', evalNs);

    refined.add(TriangleMesh.Create(objectToWorld, worldToObject,
                                         reverseOrientation, paramSet));
  }

  static Point NurbsEvaluateSurface(int uOrder, List<double> uKnot, int ucp,
                                    double u, int vOrder, List<double> vKnot,
                                    int vcp, double v, List<double> cp,
                                    [Vector dPdu, Vector dPdv]) {
    Float32List iso = new Float32List(Math.max(uOrder, vOrder) * 4);
    int uOffset = KnotOffset(uKnot, uOrder, ucp, u);
    int uFirstCp = uOffset - uOrder + 1;
    assert(uFirstCp >= 0 && uFirstCp + uOrder - 1 < ucp);

    for (int i = 0, j = 0; i < uOrder; ++i) {
      List pt = NurbsEvaluate(vOrder, vKnot, cp, (uFirstCp + i) * 4,
                              vcp, ucp, v);
      iso[j++] = pt[0];
      iso[j++] = pt[1];
      iso[j++] = pt[2];
      iso[j++] = pt[3];
    }

    int vOffset = KnotOffset(vKnot, vOrder, vcp, v);
    int vFirstCp = vOffset - vOrder + 1;
    assert(vFirstCp >= 0 && vFirstCp + vOrder - 1 < vcp);

    List<double> P = NurbsEvaluate(uOrder, uKnot, iso, (-uFirstCp) * 4, ucp, 1,
                                   u, dPdu);

    if (dPdv != null) {
      for (int i = 0, j = 0; i < vOrder; ++i) {
        List pt = NurbsEvaluate(uOrder, uKnot, cp, ((vFirstCp + i) * ucp) * 4,
                                ucp, 1, u);
        iso[j++] = pt[0];
        iso[j++] = pt[1];
        iso[j++] = pt[2];
        iso[j++] = pt[3];
      }

      NurbsEvaluate(vOrder, vKnot, iso, -vFirstCp * 4, vcp, 1, v, dPdv);
    }

    return new Point(P[0] / P[3], P[1] / P[3], P[2] / P[3]);
  }

  static List<double> NurbsEvaluate(int order, List<double> knot,
                                    List<double> cp, int cpi,
                                    int np, int cpStride,
                                    double t, [Vector deriv]) {
    int knotOffset = KnotOffset(knot, order, np, t);

    int cpOffset = knotOffset - order + 1;
    assert(cpOffset >= 0 && cpOffset < np);

    Float32List cpWork = new Float32List(4 * order);
    for (int i = 0; i < cpWork.length; ) {
      int j = cpi + ((cpOffset + i) * cpStride);
      cpWork[i++] = cp[j++];
      cpWork[i++] = cp[j++];
      cpWork[i++] = cp[j++];
      cpWork[i++] = cp[j++];
    }

    for (int i = 0; i < order - 2; ++i) {
      for (int j = 0, k = 0, l = 4; j < order - 1 - i; ++j, k += 4, l += 4) {
        double alpha = (knot[knotOffset + 1 + j] - t) /
                       (knot[knotOffset + 1 + j] -
                        knot[knotOffset + j + 2 - order + i]);
        assert(alpha >= 0.0 && alpha <= 1.0);

        cpWork[k] = cpWork[j] * alpha + cpWork[l] * (1 - alpha);
        cpWork[k + 1] = cpWork[j + 1] * alpha + cpWork[l + 1] * (1 - alpha);
        cpWork[k + 2] = cpWork[j + 2] * alpha + cpWork[l + 2] * (1 - alpha);
        cpWork[k + 3] = cpWork[j + 3] * alpha + cpWork[l + 3] * (1 - alpha);
      }
    }

    double alpha = (knot[knotOffset + 1] - t) /
                   (knot[knotOffset + 1] - knot[knotOffset + 0]);
    assert(alpha >= 0.0 && alpha <= 1.0);

    double vx = cpWork[0] * alpha + cpWork[4] * (1 - alpha);
    double vy = cpWork[1] * alpha + cpWork[5] * (1 - alpha);
    double vz = cpWork[2] * alpha + cpWork[6] * (1 - alpha);
    double vw = cpWork[3] * alpha + cpWork[7] * (1 - alpha);

    if (deriv != null) {
      double factor = (order - 1) / (knot[knotOffset + 1] - knot[knotOffset]);
      double dx = (cpWork[4] - cpWork[0]) * factor;
      double dy = (cpWork[5] - cpWork[1]) * factor;
      double dz = (cpWork[6] - cpWork[2]) * factor;
      double dw = (cpWork[7] - cpWork[3]) * factor;

      deriv.x = dx / vw - (vx * dw / (vw * vw));
      deriv.y = dy / vw - (vy * dw / (vw * vw));
      deriv.z = dz / vw - (vz * dw / (vw * vw));
    }

    return [vx, vy, vz, vw];
  }

  static int KnotOffset(List<double> knot, int order, int np, double t) {
    int firstKnot = order - 1;

    int knotOffset = firstKnot;
    while (t > knot[knotOffset + 1]) {
      ++knotOffset;
    }
    assert(knotOffset < np); // np == lastKnot
    assert(t >= knot[knotOffset] && t <= knot[knotOffset + 1]);
    return knotOffset;
  }

  int nu;
  int uorder;
  int nv;
  int vorder;
  double umin;
  double umax;
  double vmin;
  double vmax;
  List<double> uknot;
  List<double> vknot;
  bool isHomogeneous;
  List<double> P;

  static Nurbs Create(Transform o2w, Transform w2o,
                      bool ReverseOrientation, ParamSet params) {
    int nu = params.findOneInt('nu', -1);
    int uorder = params.findOneInt('uorder', -1);
    List<double> uknots = params.findFloat('uknots');
    assert(nu != -1 && uorder != -1 && uknots != null);
    assert(uknots.length == nu + uorder);
    double u0 = params.findOneFloat('u0', uknots[uorder - 1]);
    double u1 = params.findOneFloat('u1', uknots[nu]);

    int nv = params.findOneInt('nv', -1);
    int vorder = params.findOneInt('vorder', -1);
    List<double> vknots = params.findFloat('vknots');
    assert(nv != -1 && vorder != -1 && vknots != null);
    assert(vknots.length == nv + vorder);
    double v0 = params.findOneFloat('v0', vknots[vorder - 1]);
    double v1 = params.findOneFloat('v1', vknots[nv]);

    bool isHomogeneous = false;
    List<double> P;
    List<Point> p = params.findPoint('P');
    int npt;
    if (p != null) {
      P = new Float32List(p.length * 3);
      for (int i = 0, j = 0; i < p.length; ++i) {
        P[j++] = p[i].x;
        P[j++] = p[i].y;
        P[j++] = p[i].z;
      }
      npt = P.length ~/ 3;
    } else {
      P = params.findFloat('Pw');
      if (P == null) {
        LogError('Must provide control points via \'P\' or \'Pw\' parameter to '
                 'NURBS shape.');
        return null;
      }
      if ((P.length % 4) != 0) {
        LogError('Number of \'Pw\' control points provided to NURBS shape '
                 'must be multiple of four');
        return null;
      }

      npt = P.length ~/ 4;
      isHomogeneous = true;
    }

    if (npt != nu * nv) {
      LogError('NURBS shape was expecting $nu*$nv=${nu*nv} control points, '
               'was given ${P.length}');
        return null;
      }

      return new Nurbs(o2w, w2o, ReverseOrientation, nu, uorder, uknots,
                            u0, u1, nv, vorder, vknots, v0, v1, P,
                            isHomogeneous);
    }
}

