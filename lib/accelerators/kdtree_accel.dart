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
part of accelerators;

/**
 * Kd-Tree ray accelerator.
 */
class KdTreeAccel extends Aggregate {
  KdTreeAccel(List<Primitive> p,
              [this.isectCost = 80, this.traversalCost = 1,
              this.emptyBonus = 0.5,
              this.maxPrims = 1, this.maxDepth = -1]) {
    LogInfo('Building Kd-Tree Acceleration Structures.');
    for (int i = 0; i < p.length; ++i) {
      p[i].fullyRefine(primitives);
    }

    // Build kd-tree for accelerator
    nextFreeNode = 0;
    nAllocedNodes = 0;
    if (maxDepth <= 0) {
      maxDepth = (8 + 1.3 * Log2(primitives.length).toInt()).round();
    }

    bounds = new BBox();

    // Compute bounds for kd-tree construction
    List<BBox> primBounds = new List<BBox>(primitives.length);
    for (int i = 0; i < primitives.length; ++i) {
      BBox b = primitives[i].worldBound();
      bounds = BBox.Union(bounds, b);
      primBounds[i] = b;
    }

    // Allocate working memory for kd-tree construction
    List<List<_BoundEdge>> edges = new List<List<_BoundEdge>>(3);
    for (int i = 0; i < 3; ++i) {
      edges[i] = new List<_BoundEdge>(2 * primitives.length);
      for (int j = 0, len = edges[i].length; j < len; ++j) {
        edges[i][j] = new _BoundEdge();
      }
    }

    Uint32List prims0 = new Uint32List(primitives.length);
    Uint32List prims1 = new Uint32List((maxDepth + 1) * primitives.length);

    // Initialize _primNums_ for kd-tree construction
    Uint32List primNums = new Uint32List(primitives.length);
    for (int i = 0; i < primitives.length; ++i) {
      primNums[i] = i;
    }

    // Start recursive construction of kd-tree
    _buildTree(0, bounds, primBounds, primNums, primitives.length,
              maxDepth, edges, prims0, prims1, 0);
  }

  BBox worldBound() {
    return bounds;
  }

  bool canIntersect() {
    return true;
  }

  bool intersect(Ray ray, Intersection isect) {
    // Compute initial parametric range of ray inside kd-tree extent
    List<double> tmin = [0.0];
    List<double> tmax = [0.0];
    if (!bounds.intersectP(ray, tmin, tmax)) {
      return false;
    }

    // Prepare to traverse kd-tree for ray
    Vector invDir = new Vector(1.0 / ray.direction.x,
                               1.0 / ray.direction.y,
                               1.0 / ray.direction.z);

    const int MAX_TODO = 64;
    List<_KdToDo> todo = new List<_KdToDo>(MAX_TODO);
    int todoPos = 0;

    // Traverse kd-tree nodes in order for ray
    bool hit = false;
    int nodeIndex = 0;
    int firstChildIndex = 0;
    int secondChildIndex = 0;
    _KdAccelNode node = nodes[nodeIndex];
    while (node != null) {
      // Bail out if we found a hit closer than the current node
      if (ray.maxDistance < tmin[0]) {
        break;
      }

      if (!node.isLeaf()) {
        // Process kd-tree interior node

        // Compute parametric distance along ray to split plane
        int axis = node.splitAxis();
        double tplane = (node.splitPos() - ray.origin[axis]) * invDir[axis];

        // Get node children pointers for ray
        _KdAccelNode firstChild;
        _KdAccelNode secondChild;
        bool belowFirst = (ray.origin[axis] <  node.splitPos()) ||
                         (ray.origin[axis] == node.splitPos() &&
                          ray.direction[axis] <= 0);

        if (belowFirst) {
          firstChildIndex = nodeIndex + 1;
          firstChild = nodes[firstChildIndex];
          secondChildIndex = node.aboveChild();
          secondChild = nodes[secondChildIndex];
        } else {
          firstChildIndex = node.aboveChild();
          firstChild = nodes[firstChildIndex];
          secondChildIndex = nodeIndex + 1;
          secondChild = nodes[secondChildIndex];
        }

        // Advance to next child node, possibly enqueue other child
        if (tplane > tmax[0] || tplane <= 0) {
          node = firstChild;
          nodeIndex = firstChildIndex;
        } else if (tplane < tmin[0]) {
          node = secondChild;
          nodeIndex = secondChildIndex;
        } else {
          // Enqueue _secondChild_ in todo list
          if (todo[todoPos] == null) {
            todo[todoPos] = new _KdToDo();
          }
          todo[todoPos].node = secondChild;
          todo[todoPos].nodeIndex = secondChildIndex;
          todo[todoPos].tmin = tplane;
          todo[todoPos].tmax = tmax[0];
          ++todoPos;
          node = firstChild;
          nodeIndex = firstChildIndex;
          tmax[0] = tplane;
        }
      } else {
        // Check for intersections inside leaf node
        int nPrimitives = node.nPrimitives();
        if (nPrimitives == 1) {
            Primitive prim = primitives[node.onePrimitive];
            // Check one primitive inside leaf node
            if (prim.intersect(ray, isect)) {
              hit = true;
            }
        } else {
          List<int> prims = node.primitives;
          for (int i = 0; i < nPrimitives; ++i) {
            Primitive prim = primitives[prims[i]];
            // Check one primitive inside leaf node
            if (prim.intersect(ray, isect)) {
              hit = true;
            }
          }
        }

        // Grab next node to process from todo list
        if (todoPos > 0) {
          --todoPos;
          node = todo[todoPos].node;
          nodeIndex = todo[todoPos].nodeIndex;
          tmin[0] = todo[todoPos].tmin;
          tmax[0] = todo[todoPos].tmax;
        } else {
          break;
        }
      }
    }

    return hit;
  }

  bool intersectP(Ray ray) {
    // Compute initial parametric range of ray inside kd-tree extent
    List<double> tmin = [0.0];
    List<double> tmax = [0.0];
    if (!bounds.intersectP(ray, tmin, tmax)) {
      return false;
    }

    // Prepare to traverse kd-tree for ray
    Vector invDir = new Vector(1.0 / ray.direction.x,
                               1.0 / ray.direction.y,
                               1.0 / ray.direction.z);

    const int MAX_TODO = 64;
    List<_KdToDo> todo = new List<_KdToDo>(MAX_TODO);
    int todoPos = 0;

    // Traverse kd-tree nodes in order for ray
    bool hit = false;
    int nodeIndex = 0;
    int firstChildIndex = 0;
    int secondChildIndex = 0;
    _KdAccelNode node = nodes[nodeIndex];
    while (node != null) {
      // Bail out if we found a hit closer than the current node
      if (ray.maxDistance < tmin[0]) {
        break;
      }

      if (!node.isLeaf()) {
        // Process kd-tree interior node

        // Compute parametric distance along ray to split plane
        int axis = node.splitAxis();
        double tplane = (node.splitPos() - ray.origin[axis]) * invDir[axis];

        // Get node children pointers for ray
        _KdAccelNode firstChild;
        _KdAccelNode secondChild;
        bool belowFirst = (ray.origin[axis] <  node.splitPos()) ||
                          (ray.origin[axis] == node.splitPos() &&
                          ray.direction[axis] <= 0);

        if (belowFirst) {
          firstChildIndex = nodeIndex + 1;
          firstChild = nodes[firstChildIndex];
          secondChildIndex = node.aboveChild();
          secondChild = nodes[secondChildIndex];
        } else {
          firstChildIndex = node.aboveChild();
          firstChild = nodes[firstChildIndex];
          secondChildIndex = nodeIndex + 1;
          secondChild = nodes[secondChildIndex];
        }

        // Advance to next child node, possibly enqueue other child
        if (tplane > tmax[0] || tplane <= 0) {
          node = firstChild;
          nodeIndex = firstChildIndex;
        } else if (tplane < tmin[0]) {
          node = secondChild;
          nodeIndex = secondChildIndex;
        } else {
          // Enqueue _secondChild_ in todo list
          if (todo[todoPos] == null) {
            todo[todoPos] = new _KdToDo();
          }
          todo[todoPos].node = secondChild;
          todo[todoPos].nodeIndex = secondChildIndex;
          todo[todoPos].tmin = tplane;
          todo[todoPos].tmax = tmax[0];
          ++todoPos;
          node = firstChild;
          nodeIndex = firstChildIndex;
          tmax[0] = tplane;
        }
      } else {
        // Check for intersections inside leaf node
        int nPrimitives = node.nPrimitives();
        if (nPrimitives == 1) {
          Primitive prim = primitives[node.onePrimitive];
          // Check one primitive inside leaf node
          if (prim.intersectP(ray)) {
            return true;
          }
        } else {
          List<int> prims = node.primitives;
          for (int i = 0; i < nPrimitives; ++i) {
            Primitive prim = primitives[prims[i]];
            // Check one primitive inside leaf node
            if (prim.intersectP(ray)) {
              return true;
            }
          }
        }

        // Grab next node to process from todo list
        if (todoPos > 0) {
          --todoPos;
          node = todo[todoPos].node;
          nodeIndex = todo[todoPos].nodeIndex;
          tmin[0] = todo[todoPos].tmin;
          tmax[0] = todo[todoPos].tmax;
        } else {
          break;
        }
      }
    }

    return false;
  }

  void _buildTree(int nodeNum, BBox nodeBounds,
      List<BBox> allPrimBounds, List<int> primNums, int nPrimitives, int depth,
      List<List<_BoundEdge>> edges, List<int> prims0, List<int> prims1,
      int prims1Offset,
      [int badRefines = 0]) {
    assert(nodeNum == nextFreeNode);

    // Get next free node from _nodes_ array
    if (nextFreeNode == nAllocedNodes) {
      int nAlloc = Math.max(2 * nAllocedNodes, 512);
      List<_KdAccelNode> n = new List<_KdAccelNode>(nAlloc);
      if (nAllocedNodes > 0) {
        for (int i = 0; i < nAllocedNodes; ++i) {
          n[i] = nodes[i];
        }

        for (int i = nAllocedNodes; i < n.length; ++i) {
          n[i] = new _KdAccelNode();
        }
      } else {
        for (int i = 0; i < n.length; ++i) {
          n[i] = new _KdAccelNode();
        }
      }

      nodes = n;
      nAllocedNodes = nAlloc;
    }
    ++nextFreeNode;

    // Initialize leaf node if termination criteria met
    if (nPrimitives <= maxPrims || depth == 0) {
      nodes[nodeNum].initLeaf(primNums, nPrimitives);
      return;
    }

    // Initialize interior node and continue recursion

    // Choose split axis position for interior node
    int bestAxis = -1;
    int bestOffset = -1;
    double bestCost = INFINITY;
    double oldCost = isectCost * nPrimitives.toDouble();
    double totalSA = nodeBounds.surfaceArea();
    double invTotalSA = 1.0 / totalSA;
    Vector d = nodeBounds.pMax - nodeBounds.pMin;

    // Choose which axis to split along
    int axis = nodeBounds.maximumExtent();
    int retries = 0;

    //retrySplit:
    while (true) {
      // Initialize edges for _axis_
      for (int i = 0; i < nPrimitives; ++i) {
        int pn = primNums[i];
        BBox bbox = allPrimBounds[pn];
        edges[axis][2 * i] = new _BoundEdge(bbox.pMin[axis], pn, true);
        edges[axis][2 * i + 1] = new _BoundEdge(bbox.pMax[axis], pn, false);
      }

      edges[axis].sort((a, b) => a < b ? -1 : 1);

      // Compute cost of all splits for _axis_ to find best
      int nBelow = 0;
      int nAbove = nPrimitives;
      for (int i = 0; i < 2 * nPrimitives; ++i) {
        if (edges[axis][i].type == _BoundEdge.END) {
          --nAbove;
        }

        double edget = edges[axis][i].t;
        if (edget > nodeBounds.pMin[axis] &&
            edget < nodeBounds.pMax[axis]) {
          // Compute cost for split at _i_th edge
          int otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;
          double belowSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                               (edget - nodeBounds.pMin[axis]) *
                               (d[otherAxis0] + d[otherAxis1]));
          double aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                               (nodeBounds.pMax[axis] - edget) *
                               (d[otherAxis0] + d[otherAxis1]));
          double pBelow = belowSA * invTotalSA;
          double pAbove = aboveSA * invTotalSA;
          double eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0.0;
          double cost = traversalCost +
                       isectCost * (1.0 - eb) *
                       (pBelow * nBelow + pAbove * nAbove);

          // Update best split if this is lowest cost so far
          if (cost < bestCost) {
            bestCost = cost;
            bestAxis = axis;
            bestOffset = i;
          }
        }

        if (edges[axis][i].type == _BoundEdge.START) {
          ++nBelow;
        }
      }

      assert(nBelow == nPrimitives && nAbove == 0);

      // Create leaf if no good splits were found
      if (bestAxis == -1 && retries < 2) {
        ++retries;
        axis = (axis + 1) % 3;
        continue;
      }

      break;
    }

    if (bestCost > oldCost) {
      ++badRefines;
    }

    if ((bestCost > 4.0 * oldCost && nPrimitives < 16) ||
        bestAxis == -1 || badRefines == 3) {
      nodes[nodeNum].initLeaf(primNums, nPrimitives);
      return;
    }

    // Classify primitives with respect to split
    int n0 = 0, n1 = 0;
    for (int i = 0; i < bestOffset; ++i) {
      if (edges[bestAxis][i].type == _BoundEdge.START) {
        prims0[n0++] = edges[bestAxis][i].primNum;
      }
    }

    for (int i = bestOffset + 1; i < 2 * nPrimitives; ++i) {
      if (edges[bestAxis][i].type == _BoundEdge.END) {
        prims1[prims1Offset + n1++] = edges[bestAxis][i].primNum;
      }
    }

    // Recursively initialize children nodes
    double tsplit = edges[bestAxis][bestOffset].t;

    BBox bounds0 = nodeBounds, bounds1 = nodeBounds;
    bounds0.pMax[bestAxis] = bounds1.pMin[bestAxis] = tsplit;
    _buildTree(nodeNum + 1, bounds0,
               allPrimBounds, prims0, n0, depth - 1, edges,
               prims0, prims1, nPrimitives, badRefines);

    int aboveChild = nextFreeNode;
    nodes[nodeNum].initInterior(bestAxis, aboveChild, tsplit);
    _buildTree(aboveChild, bounds1, allPrimBounds, prims1, n1,
               depth - 1, edges, prims0, prims1, nPrimitives, badRefines);
  }

  static KdTreeAccel Create(List<Primitive> prims, ParamSet ps) {
    int isectCost = ps.findOneInt('intersectcost', 80);
    int travCost = ps.findOneInt('traversalcost', 1);
    double emptyBonus = ps.findOneFloat('emptybonus', 0.5);
    int maxPrims = ps.findOneInt('maxprims', 1);
    int maxDepth = ps.findOneInt('maxdepth', -1);
    return new KdTreeAccel(prims, isectCost, travCost,
        emptyBonus, maxPrims, maxDepth);
  }

  int isectCost;
  int traversalCost;
  int maxPrims;
  int maxDepth;
  double emptyBonus;
  List<Primitive> primitives = [];
  List<_KdAccelNode> nodes = [];
  int nAllocedNodes;
  int nextFreeNode;
  BBox bounds;
}

class _KdToDo {
  _KdAccelNode node;
  int nodeIndex;
  double tmin;
  double tmax;
}

class _KdAccelNode {
  void initLeaf(List<int> primNums, int np) {
    flags = 3;
    flags |= (np << 2);
    // Store primitive ids for leaf node
    if (np == 0) {
      onePrimitive = 0;
    } else if (np == 1) {
      onePrimitive = primNums[0];
    } else {
      primitives = new Uint32List(np);
      for (int i = 0; i < np; ++i) {
        primitives[i] = primNums[i];
      }
    }
  }

  void initInterior(int axis, int ac, double s) {
    split = s;
    flags = axis;
    flags |= (ac << 2);
  }

  double splitPos() {
    return split;
  }

  int nPrimitives() {
    return flags >> 2;
  }

  int splitAxis() {
    return flags & 3;
  }

  bool isLeaf() {
    return (flags & 3) == 3;
  }

  int aboveChild() {
    return flags >> 2;
  }

  double split;
  int onePrimitive;
  Uint32List primitives;

  int flags;
}


class _BoundEdge {
  _BoundEdge([double tt = 0.0, int pn = 0, bool starting = false]) {
    t = tt;
    primNum = pn;
    type = starting ? START : END;
  }

  bool operator <(_BoundEdge e) {
    if (t == e.t) {
      return type < e.type;
    }
    else return t < e.t;
  }

  double t;
  int primNum;
  int type;

  static const int START = 1;
  static const int END = 0;
}
