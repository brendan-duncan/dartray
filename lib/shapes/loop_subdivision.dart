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

class LoopSubdivision extends Shape {
  LoopSubdivision(Transform o2w, Transform w2o, bool ro,
             int nfaces, int nvertices, List<int> vertexIndices,
             List<Point> P, this.nLevels) :
    super(o2w, w2o, ro),
    vertices = new List<_SDVertex>(nvertices),
    faces = new List<_SDFace>(nfaces) {

    for (int i = 0; i < nvertices; ++i) {
      vertices[i] = new _SDVertex(P[i]);
    }

    // Set face to vertex pointers
    for (int i = 0, j = 0; i < nfaces; ++i) {
      _SDFace f = new _SDFace();
      faces[i] = f;

      _SDVertex v = vertices[vertexIndices[j++]];
      f.v[0] = v;
      v.startFace = f;

      v = vertices[vertexIndices[j++]];
      f.v[1] = v;
      v.startFace = f;

      v = vertices[vertexIndices[j++]];
      f.v[2] = v;
      v.startFace = f;
    }

    _SDEdgeMap edges = new _SDEdgeMap();
    for (int i = 0; i < nfaces; ++i) {
      _SDFace f = faces[i];

      for (int ei = 0; ei < 3; ++ei) {
        _SDVertex v0 = f.v[ei];
        _SDVertex v1 = f.v[(ei + 1) % 3];
        _SDEdge edge = edges.getEdge(v0, v1);

        if (edge == null) {
          edge = new _SDEdge(v0, v1);
          edge.f[0] = f;
          edge.f0edgeNum = ei;
          edges.setEdge(v0, v1, edge);
        } else {
          edge.f[0].f[edge.f0edgeNum] = f;
          f.f[ei] = edge.f[0];
        }
      }
    }

    // Finish vertex initialization
    for (int i = 0; i < nvertices; ++i) {
      _SDVertex v = vertices[i];
      _SDFace f = v.startFace;

      do {
        f = f.nextFace(v);
      } while (f != null && f != v.startFace);

      v.boundary = (f == null);

      int val = v.valence();
      if (!v.boundary && val == 6) {
        v.regular = true;
      } else if (v.boundary && val == 4) {
        v.regular = true;
      } else {
        v.regular = false;
      }
    }
  }

  bool canIntersect() {
    return false;
  }

  void refine(List<Shape> refined) {
    List<_SDFace> f = faces;
    List<_SDVertex> v = vertices;

    for (int i = 0; i < nLevels; ++i) {
      // Update _f_ and _v_ for next level of subdivision
      List<_SDFace> newFaces = [];
      List<_SDVertex> newVertices = [];

      // Allocate next level of children in mesh tree
      for (int j = 0; j < v.length; ++j) {
        v[j].child = new _SDVertex(null);
        v[j].child.regular = v[j].regular;
        v[j].child.boundary = v[j].boundary;
        newVertices.add(v[j].child);
      }

      for (int j = 0; j < f.length; ++j) {
        for (int k = 0; k < 4; ++k) {
          f[j].children[k] = new _SDFace();
          newFaces.add(f[j].children[k]);
        }
      }

      // Update vertex positions and create new edge vertices

      // Update vertex positions for even vertices
      for (int j = 0; j < v.length; ++j) {
        if (!v[j].boundary) {
          // Apply one-ring rule for even vertex
          if (v[j].regular) {
            v[j].child.P = weightOneRing(v[j], 1.0 / 16.0);
          } else {
            v[j].child.P = weightOneRing(v[j], beta(v[j].valence()));
          }
        } else {
          // Apply boundary rule for even vertex
          v[j].child.P = weightBoundary(v[j], 1.0 / 8.0);
        }
      }

      // Compute new odd edge vertices
      _SDEdgeMap edgeVerts = new _SDEdgeMap();

      for (int j = 0; j < f.length; ++j) {
        _SDFace face = f[j];
        for (int k = 0; k < 3; ++k) {
          // Compute odd vertex on _k_th edge
          _SDVertex fv1 = face.v[k];
          _SDVertex fv2 = face.v[(k + 1) % 3];
          _SDVertex vert = edgeVerts.getEdge(fv1, fv2);

          if (vert == null) {
            // Create and initialize new odd vertex
            vert = new _SDVertex(null);
            newVertices.add(vert);

            vert.regular = true;
            vert.boundary = (face.f[k] == null);
            vert.startFace = face.children[3];

            // Apply edge rules to compute new vertex position
            if (vert.boundary) {
              vert.P = fv1.P * 0.5 + fv2.P * 0.5;
            } else {
              vert.P = fv1.P * (3.0 / 8.0) + fv2.P * (3.0 / 8.0);
              vert.P += face.otherVert(fv1, fv2).P * (1.0 / 8.0);
              vert.P += face.f[k].otherVert(fv1, fv2).P * (1.0 / 8.0);
            }

            edgeVerts.setEdge(fv1, fv2, vert);
          }
        }
      }

      // Update new mesh topology

      // Update even vertex face pointers
      for (int j = 0; j < v.length; ++j) {
        _SDVertex vert = v[j];
        int vertNum = vert.startFace.vnum(vert);
        vert.child.startFace = vert.startFace.children[vertNum];
      }

      // Update face neighbor pointers
      for (int j = 0; j < f.length; ++j) {
        _SDFace face = f[j];
        for (int k = 0; k < 3; ++k) {
          // Update children _f_ pointers for siblings
          face.children[3].f[k] = face.children[(k + 1) % 3];
          face.children[k].f[(k + 1) % 3] = face.children[3];

          // Update children _f_ pointers for neighbor children
          _SDFace f2 = face.f[k];
          face.children[k].f[k] = f2 != null ?
                                  f2.children[f2.vnum(face.v[k])] :
                                  null;

          f2 = face.f[(k + 2) % 3];
          face.children[k].f[(k + 2) % 3] = f2 != null ?
                                            f2.children[f2.vnum(face.v[k])] :
                                            null;
        }
      }

      // Update face vertex pointers
      for (int j = 0; j < f.length; ++j) {
        _SDFace face = f[j];
        for (int k = 0; k < 3; ++k) {
          // Update child vertex pointer to new even vertex
          face.children[k].v[k] = face.v[k].child;

          // Update child vertex pointer to new odd vertex
          _SDVertex fv1 = face.v[k];
          _SDVertex fv2 = face.v[(k + 1) % 3];
          _SDVertex vert = edgeVerts.getEdge(fv1, fv2);

          face.children[k].v[(k + 1) % 3] = vert;
          face.children[(k + 1) % 3].v[k] = vert;
          face.children[3].v[k] = vert;
        }
      }

      // Prepare for next level of subdivision
      f = newFaces;
      v = newVertices;
    }

    // Push vertices to limit surface
    List<Point> Plimit = new List<Point>(v.length);
    for (int i = 0; i < v.length; ++i) {
      if (v[i].boundary) {
        Plimit[i] = weightBoundary(v[i], 1.0 / 5.0);
      } else {
        Plimit[i] = weightOneRing(v[i], gamma(v[i].valence()));
      }
    }

    for (int i = 0; i < v.length; ++i) {
      v[i].P = Plimit[i];
    }

    // Compute vertex tangents on limit surface
    List<Normal> Ns = [];
    List<Point> Pring = new List<Point>();
    Pring.length = 16;

    for (int i = 0; i < v.length; ++i) {
      _SDVertex vert = v[i];
      Vector S = new Vector();
      Vector T = new Vector();
      int valence = vert.valence();
      if (valence > Pring.length) {
        Pring.length = valence;
      }
      vert.oneRing(Pring, 0);
      if (!vert.boundary) {
        // Compute tangents of interior face
        for (int k = 0; k < valence; ++k) {
          S += Pring[k] * (Math.cos(2.0 * Math.PI * k / valence));
          T += Pring[k] * (Math.sin(2.0 * Math.PI * k / valence));
        }
      } else {
        // Compute tangents of boundary face
        S = Pring[valence - 1] - Pring[0];
        if (valence == 2) {
          T = Pring[0] + Pring[1] - vert.P * 2.0;
        } else if (valence == 3) {
          T = Pring[1] - vert.P;
        } else if (valence == 4) { // regular
          T = Pring[0] * -1.0 + Pring[1] * 2.0 + Pring[2] * 2.0 +
              Pring[3] * -1.0 + vert.P * -2.0;
        } else {
          double theta = Math.PI / (valence - 1);
          T = (Pring[0] + Pring[valence - 1]) * Math.sin(theta);
          for (int k = 1; k < valence - 1; ++k) {
            double wt = (2 * Math.cos(theta) - 2) * Math.sin((k) * theta);
            T += Pring[k] * wt;
          }
          T = -T;
        }
      }

      Ns.add(new Normal.from(Vector.Cross(S, T)));
    }

    // Create _TriangleMesh_ from subdivision mesh
    int ntris = f.length;
    Uint32List verts = new Uint32List(3 * ntris);
    int vi = 0;
    int totVerts = v.length;
    Map<_SDVertex, int> usedVerts = {};
    for (int i = 0; i < totVerts; ++i) {
      usedVerts[v[i]] = i;
    }

    for (int i = 0; i < ntris; ++i) {
      for (int j = 0; j < 3; ++j) {
        verts[vi++] = usedVerts[f[i].v[j]];
      }
    }

    ParamSet paramSet = new ParamSet();
    paramSet.addInt('indices', verts);
    paramSet.addPoint('P', Plimit);
    paramSet.addNormal('N', Ns);

    refined.add(TriangleMesh.Create(objectToWorld,
                worldToObject, reverseOrientation, paramSet));
  }

  BBox objectBound() {
    BBox b = new BBox();
    for (int i = 0; i < vertices.length; i++) {
      b = BBox.UnionPoint(b, vertices[i].P);
    }
    return b;
  }

  BBox worldBound() {
    BBox b = new BBox();
    for (int i = 0; i < vertices.length; i++) {
      b = BBox.UnionPoint(b, objectToWorld.transformPoint(vertices[i].P));
    }
    return b;
  }

  static double beta(int valence) {
    if (valence == 3) {
      return 3.0 / 16.0;
    }
    return 3.0 / (8.0 * valence);
  }

  static Point weightOneRing(_SDVertex vert, double beta) {
    // Put _vert_ one-ring in _Pring_
    int valence = vert.valence();
    List<Point> Pring = new List<Point>(valence);
    vert.oneRing(Pring);
    Point P = vert.P * (1.0 - valence * beta);
    for (int i = 0; i < valence; ++i) {
      P += Pring[i] * beta;
    }
    return P;
  }

  static Point weightBoundary(_SDVertex vert, double beta) {
    // Put _vert_ one-ring in _Pring_
    int valence = vert.valence();
    List<Point> Pring = new List<Point>(valence);
    vert.oneRing(Pring);
    Point P = vert.P * (1.0 - 2.0 * beta);
    P += Pring[0] * beta;
    P += Pring[valence - 1] * beta;
    return P;
  }

  static double gamma(int valence) {
    return 1.0 / (valence + 3.0 / (8.0 * beta(valence)));
  }

  int nLevels;
  final List<_SDVertex> vertices;
  final List<_SDFace> faces;

  static LoopSubdivision Create(Transform o2w, Transform w2o,
                                bool reverseOrientation, ParamSet params) {
      int nlevels = params.findOneInt('nlevels', 1);
      List<int> vi = params.findInt('indices');
      List<Point> P = params.findPoint('P');
      if (vi == null || P == null) {
        return null;
      }

      return new LoopSubdivision(o2w, w2o, reverseOrientation,
                                 vi.length ~/ 3, P.length, vi,
                                 P, nlevels);
  }
}

class _SDEdgeMap {
  getEdge(a, b) {
    if (_edgeMap.containsKey(a)) {
      if (_edgeMap[a].containsKey(b)) {
        return _edgeMap[a][b];
      }
    }
    if (_edgeMap.containsKey(b)) {
      if (_edgeMap[b].containsKey(a)) {
        return _edgeMap[b][a];
      }
    }
    return null;
  }

  void setEdge(a, b, v) {
    if (!_edgeMap.containsKey(a)) {
      _edgeMap[a] = new Map();
    }
    _edgeMap[a][b] = v;
  }

  Map _edgeMap = {};
}

class _SDVertex {
  _SDVertex(this.P) :
    regular = false,
    boundary = false;

  int valence() {
    _SDFace f = startFace;
    if (!boundary) {
      // Compute valence of interior vertex
      int nf = 1;
      while ((f = f.nextFace(this)) != startFace) {
        ++nf;
      }
      return nf;
    } else {
      // Compute valence of boundary vertex
      int nf = 1;
      while ((f = f.nextFace(this)) != null) {
        ++nf;
      }
      f = startFace;
      while ((f = f.prevFace(this)) != null) {
        ++nf;
      }
      return nf + 1;
    }
  }

  void oneRing(List<Point> p, [int pi = 0]) {
    if (!boundary) {
      // Get one-ring vertices for interior vertex
      _SDFace face = startFace;
      do {
        p[pi++] = face.nextVert(this).P;
        face = face.nextFace(this);
      } while (face != startFace);
    } else {
      // Get one-ring vertices for boundary vertex
      _SDFace face = startFace;
      _SDFace f2;
      while ((f2 = face.nextFace(this)) != null) {
        face = f2;
      }
      p[pi++] = face.nextVert(this).P;
      do {
        p[pi++] = face.prevVert(this).P;
        face = face.prevFace(this);
      } while (face != null);
    }
  }

  Point P;
  _SDFace startFace;
  _SDVertex child;
  bool regular;
  bool boundary;
}

class _SDFace {
  int vnum(_SDVertex vert) {
    for (int i = 0; i < 3; ++i) {
      if (v[i] == vert) {
        return i;
      }
    }
    LogSevere('Basic logic error in SDFace::vnum()');
    return -1;
  }

  _SDFace nextFace(_SDVertex vert) {
    return f[vnum(vert)];
  }

  _SDFace prevFace(_SDVertex vert) {
    return f[(vnum(vert) + 2) % 3];
  }

  _SDVertex nextVert(_SDVertex vert) {
    return v[(vnum(vert) + 1) % 3];
  }

  _SDVertex prevVert(_SDVertex vert) {
    return v[(vnum(vert) + 2) % 3];
  }

  _SDVertex otherVert(_SDVertex v0, _SDVertex v1) {
    for (int i = 0; i < 3; ++i) {
      if (v[i] != v0 && v[i] != v1) {
        return v[i];
      }
    }
    LogSevere('Basic logic error in SDVertex::otherVert()');
    return null;
  }

  final List<_SDVertex> v = [null, null, null];
  final List<_SDFace> f = [null, null, null];
  final List<_SDFace> children = [null, null, null, null];
}

class _SDEdge {
  _SDEdge(_SDVertex v0, _SDVertex v1) {
    v[0] = v0;
    v[1] = v1;
    f0edgeNum = -1;
  }

  final List<_SDVertex> v = [null, null];
  final List<_SDFace> f = [null, null];
  int f0edgeNum;
}

