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

class Heightfield extends Shape {
  Heightfield(Transform o2w, Transform w2o, bool ro, this.nx, this.ny,
              this.z) :
    super(o2w, w2o, ro);

  bool canIntersect() {
    return false;
  }

  void refine(List<Shape> refined) {
    int ntris = 2 * (nx - 1) * (ny - 1);
    Uint32List verts = new Uint32List(3 * ntris);
    List<Point> P = new List<Point>(nx * ny);
    Float32List uvs = new Float32List(2 * nx * ny);
    final int nverts = nx * ny;

    // Compute heightfield vertex positions
    for (int y = 0, pi = 0, ui = 0; y < ny; ++y) {
      for (int x = 0; x < nx; ++x, ++pi, ui += 2) {
        uvs[ui] = x / (nx - 1);
        uvs[ui + 1] = y / (ny - 1);
        P[pi] = new Point(uvs[ui], uvs[ui + 1], z[pi]);
      }
    }

    // Fill in heightfield vertex offset array
    int vp = 0;
    int VERT(int x, int y) => (x + y * nx);
    for (int y = 0; y < ny - 1; ++y) {
      for (int x = 0; x < nx - 1; ++x) {
        verts[vp++] = VERT(x, y);
        verts[vp++] = VERT(x + 1, y);
        verts[vp++] = VERT(x + 1, y + 1);

        verts[vp++] = VERT(x, y);
        verts[vp++] = VERT(x + 1, y + 1);
        verts[vp++] = VERT(x, y + 1);
      }
    }

    ParamSet paramSet = new ParamSet();
    paramSet.addInt('indices', verts);
    paramSet.addFloat('uv', uvs);
    paramSet.addPoint('P', P);
    refined.add(TriangleMesh.Create(objectToWorld, worldToObject,
                                         reverseOrientation, paramSet));
  }

  BBox objectBound() {
    double minz = z[0];
    double maxz = z[0];
    for (int i = 1; i < nx * ny; ++i) {
      if (z[i] < minz) {
        minz = z[i];
      }
      if (z[i] > maxz) {
        maxz = z[i];
      }
    }

    return new BBox(new Point(0.0, 0.0, minz), new Point(1.0, 1.0, maxz));
  }

  List<double> z;
  int nx;
  int ny;

  static Heightfield Create(Transform o2w, Transform w2o,
                            bool reverseOrientation, ParamSet params) {
    int nu = params.findOneInt('nu', -1);
    int nv = params.findOneInt('nv', -1);
    List<double> Pz = params.findFloat('Pz');
    assert(Pz.length == nu * nv);
    assert(nu != -1 && nv != -1 && Pz != null);
    return new Heightfield(o2w, w2o, reverseOrientation, nu, nv, Pz);
  }
}
