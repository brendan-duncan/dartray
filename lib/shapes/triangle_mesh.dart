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

class TriangleMesh extends Shape {
  TriangleMesh(Transform o2w, Transform w2o, bool ro,
               this.ntris, this.nverts, this.vertexIndex,
               List<Point> P, this.n, this.s,
               this.uvs, this.alphaTexture) :
    super(o2w, w2o, ro) {
    _p = new Float32List(nverts * 3);
    // Transform mesh vertices to world space
    for (int i = 0, j = 0; i < nverts; ++i) {
      Point p = objectToWorld.transformPoint(P[i]);
      _p[j++] = p.x;
      _p[j++] = p.y;
      _p[j++] = p.z;
    }
  }

  Point point(int index) {
    int i3 = index * 3;
    return new Point(_p[i3], _p[i3 + 1], _p[i3 + 2]);
  }

  List<Point> _triangle = [new Point(), new Point(), new Point()];
  List<Point> triangle(int v1, int v2, int v3) {
    int i3 = v1 * 3;
    _triangle[0].data[0] = _p[i3];
    _triangle[0].data[1] = _p[i3 + 1];
    _triangle[0].data[2] = _p[i3 + 2];
    i3 = v2 * 3;
    _triangle[1].data[0] = _p[i3];
    _triangle[1].data[1] = _p[i3 + 1];
    _triangle[1].data[2] = _p[i3 + 2];
    i3 = v3 * 3;
    _triangle[2].data[0] = _p[i3];
    _triangle[2].data[1] = _p[i3 + 1];
    _triangle[2].data[2] = _p[i3 + 2];
    return _triangle;
  }

  BBox objectBound() {
    BBox objectBounds = new BBox();
    for (int i = 0; i < nverts; i++) {
      objectBounds = BBox.UnionPoint(objectBounds,
                                     worldToObject.transformPoint(point(i)));
    }
    return objectBounds;
  }

  BBox worldBound() {
    BBox worldBounds = new BBox();
    for (int i = 0; i < nverts; i++) {
      worldBounds = BBox.UnionPoint(worldBounds, point(i));
    }
    return worldBounds;
  }

  bool canIntersect() {
    return false;
  }

  void refine(List<Shape> refined) {
    for (int i = 0; i < ntris; ++i) {
      refined.add(new Triangle(objectToWorld,
                               worldToObject, reverseOrientation,
                               this, i));
    }
  }

  int ntris;
  int nverts;
  List<int> vertexIndex;
  Float32List _p;
  List<Normal> n;
  List<Vector> s;
  List<double> uvs;
  Texture alphaTexture;


  static TriangleMesh Create(Transform o2w, Transform w2o,
                             bool reverseOrientation, ParamSet params,
                             [Map<String, Texture> floatTextures = null]) {
    List<int> vi = params.findInt('indices');
    List<Point> P = params.findPoint('P');
    List<double> uvs = params.findFloat('uv');
    if (uvs == null) {
      uvs = params.findFloat('st');
    }

    if (vi == null || P == null) {
      return null;
    }

    // Replace the List<*> with a more compact typed data lists.
    vi = new Uint32List.fromList(vi);

    if (uvs != null) {
      uvs = new Float32List.fromList(uvs);
    }

    bool discardDegnerateUVs = params.findOneBool('discarddegenerateUVs', false);

    if (uvs != null) {
      int npi = P.length;
      int nuvi = uvs.length;

      if (uvs.length < 2 * P.length) {
        LogWarning('Not enough of \'uv\'s for triangle mesh. '
                   'Expencted ${2*npi}, found ${nuvi}.  Discarding.');
        uvs = null;
      } else if (uvs.length > 2 * P.length) {
        LogWarning('More \'uv\'s provided than will be used for triangle '
                   'mesh.  (${2*npi} expcted, ${nuvi} found)');
      }
    }


    List<Vector> S = params.findVector('S');
    if (S != null && S.length != P.length) {
      LogWarning('Number of \'S\'s for triangle mesh must match \'P\'s');
      S = null;
    }

    List<Normal> N = params.findNormal('N');
    if (N != null && N.length != P.length) {
      LogWarning('Number of \'N\'s for triangle mesh must match \'P\'s');
      N = null;
    }

    if (discardDegnerateUVs && uvs != null && N != null) {
      int nvi = N.length;
      // if there are normals, check for bad uv's that
      // give degenerate mappings; discard them if so
      int vp = 0;
      for (int i = 0; i < nvi; i += 3, vp += 3) {
        double area = 0.5 * Vector.Cross(P[vi[vp + 0]] - P[vi[vp + 1]],
                                        P[vi[vp + 2]] - P[vi[vp + 1]]).length();
        if (area < 1.0e-7) {
          continue; // ignore degenerate tris.
        }

        if ((uvs[2 * vi[vp + 0]] == uvs[2 * vi[vp + 1]] &&
             uvs[2 * vi[vp + 0] + 1] == uvs[2 * vi[vp + 1] + 1]) ||
            (uvs[2 * vi[vp + 1]] == uvs[2 * vi[vp + 2]] &&
             uvs[2 * vi[vp + 1] + 1] == uvs[2 * vi[vp + 2] + 1]) ||
            (uvs[2 * vi[vp + 2]] == uvs[2 * vi[vp + 0]] &&
             uvs[2 * vi[vp + 2] + 1] == uvs[2 * vi[vp + 0] + 1])) {
          LogWarning('Degenerate uv coordinates in triangle mesh. '
                     'Discarding all uvs.');
          uvs = null;
          break;
        }
      }
    }

    int npi = P.length;
    for (int i = 0, nvi = vi.length; i < nvi; ++i) {
      if (vi[i] >= npi) {
        LogWarning('TriangleMesh has out of-bounds vertex index '
                   '${vi[i]} ($npi \'P\' values were given');
        return null;
      }
    }

    Texture alphaTex = null;
    String alphaTexName = params.findTexture('alpha');

    if (alphaTexName != '') {
      if (floatTextures.containsKey(alphaTexName)) {
        alphaTex = floatTextures[alphaTexName];
      } else {
        LogWarning('Couldn\'t find float texture \'${alphaTexName}\' '
                   'for \'alpha\' parameter');
      }
    } else if (params.findOneFloat('alpha', 1.0) == 0.0) {
      alphaTex = new ConstantTexture(0.0);
    }

    return new TriangleMesh(o2w, w2o, reverseOrientation,
                            vi.length ~/ 3, P.length, vi, P, N, S, uvs,
                            alphaTex);
  }
}
