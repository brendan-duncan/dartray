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
 *
 * Accelerator "kdtree" <Parameters>
 * ----------------------------------------------------------------------------
 *  **Type** | **Name**      | **Default** | **Description**
 * ----------------------------------------------------------------------------
 *  integer  | intersectcost | 80          |
 * ----------------------------------------------------------------------------
 *  integer  | traversalcost | 1           |
 * ----------------------------------------------------------------------------
 *  float    | emptybonus    | 0.5         |
 * ----------------------------------------------------------------------------
 *  integer  | maxprims      | 1           |
 * ----------------------------------------------------------------------------
 *  integer  | maxdepth      | -1          |
 * ----------------------------------------------------------------------------
 */
class KdTreeAccel extends Aggregate {
  KdTreeAccel(List<Primitive> p, [this.isectCost = 80, this.traversalCost = 1,
              this.emptyBonus = 0.5, this.maxPrims = 1, this.maxDepth = -1]) {
    Stats.KDTREE_STARTED_CONSTRUCTION(this, p.length);
    LogInfo('Building Kd-Tree Acceleration Structures.');
    for (int i = 0; i < p.length; ++i) {
      p[i].fullyRefine(primitives);
    }

    // Build kd-tree for accelerator
    nextFreeNode = 0;
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
    List<List<_BoundEdge>> edges = new List(3);
    for (int i = 0; i < 3; ++i) {
      edges[i] = new List<_BoundEdge>(2 * primitives.length);
      for (int j = 0, len = edges[i].length; j < len; ++j) {
        edges[i][j] = new _BoundEdge();
      }
    }

    Uint32List prims0 = new Uint32List(primitives.length);
    Uint32List prims1 = new Uint32List((maxDepth + 1) * primitives.length);
    Uint32List primNums = new Uint32List(primitives.length);

    // Initialize primNums for kd-tree construction
    for (int i = 0; i < primitives.length; ++i) {
      primNums[i] = i;
    }

    // Start recursive construction of kd-tree
    _buildTree(0, bounds, primBounds, primNums, maxDepth, edges, prims0,
               prims1);

    Stats.KDTREE_FINISHED_CONSTRUCTION(this);
  }

  BBox worldBound() {
    return bounds;
  }

  bool canIntersect() {
    return true;
  }

  bool intersect(Ray ray, Intersection isect) {
    Stats.KDTREE_INTERSECTION_TEST(this, ray);
    // Compute initial parametric range of ray inside kd-tree extent
    List<double> tmin = [0.0];
    List<double> tmax = [0.0];
    if (!bounds.intersectP(ray, tmin, tmax)) {
      Stats.KDTREE_RAY_MISSED_BOUNDS();
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
    _KdAccelNode node = nodes[nodeIndex];

    while (node != null) {
      // Bail out if we found a hit closer than the current node
      if (ray.maxDistance < tmin[0]) {
        break;
      }

      if (!node.isLeaf) {
        // Process kd-tree interior node
        Stats.KDTREE_INTERSECTION_TRAVERSED_INTERIOR_NODE(node);

        // Compute parametric distance along ray to split plane
        int axis = node.splitAxis;
        double tplane = (node.splitPos - ray.origin[axis]) * invDir[axis];

        // Get node children pointers for ray
        int firstChildIndex;
        _KdAccelNode firstChild;
        int secondChildIndex;
        _KdAccelNode secondChild;

        bool belowFirst = (ray.origin[axis] < node.splitPos) ||
                          (ray.origin[axis] == node.splitPos &&
                           ray.direction[axis] <= 0);

        if (belowFirst) {
          firstChildIndex = nodeIndex + 1;
          firstChild = nodes[firstChildIndex];

          secondChildIndex = node.aboveChild;
          secondChild = nodes[secondChildIndex];
        } else {
          firstChildIndex = node.aboveChild;
          firstChild = nodes[firstChildIndex];

          secondChildIndex = nodeIndex + 1;
          secondChild = nodes[secondChildIndex];
        }

        // Advance to next child node, possibly enqueue other child
        if (tplane > tmax[0] || tplane <= 0.0) {
          nodeIndex = firstChildIndex;
          node = firstChild;
        } else if (tplane < tmin[0]) {
          nodeIndex = secondChildIndex;
          node = secondChild;
        } else {
          // Enqueue secondChild in todo list
          if (todo[todoPos] == null) {
            todo[todoPos] = new _KdToDo();
          }

          todo[todoPos].nodeIndex = secondChildIndex;
          todo[todoPos].node = secondChild;
          todo[todoPos].tmin = tplane;
          todo[todoPos].tmax = tmax[0];
          ++todoPos;

          nodeIndex = firstChildIndex;
          node = firstChild;

          tmax[0] = tplane;
        }
      } else {
        Stats.KDTREE_INTERSECTION_TRAVERSED_LEAF_NODE(node, node.nPrimitives);
        // Check for intersections inside leaf node
        int nPrimitives = node.nPrimitives;

        if (nPrimitives == 1) {
          Primitive prim = primitives[node.primitives];
          Stats.KDTREE_INTERSECTION_PRIMITIVE_TEST(prim);
          // Check one primitive inside leaf node
          if (prim.intersect(ray, isect)) {
            Stats.KDTREE_INTERSECTION_HIT(prim);
            hit = true;
          }
        } else if (nPrimitives > 1) {
          List<int> prims = node.primitives;
          for (int i = 0; i < nPrimitives; ++i) {
            Primitive prim = primitives[prims[i]];
            // Check one primitive inside leaf node
            Stats.KDTREE_INTERSECTION_PRIMITIVE_TEST(prim);
            if (prim.intersect(ray, isect)) {
              Stats.KDTREE_INTERSECTION_HIT(prim);
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

    Stats.KDTREE_INTERSECTION_FINISHED();
    return hit;
  }

  bool intersectP(Ray ray) {
    // Compute initial parametric range of ray inside kd-tree extent
    Stats.KDTREE_INTERSECTIONP_TEST(this, ray);

    List<double> tmin = [0.0];
    List<double> tmax = [0.0];
    if (!bounds.intersectP(ray, tmin, tmax)) {
      Stats.KDTREE_RAY_MISSED_BOUNDS();
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
    _KdAccelNode node = nodes[nodeIndex];

    while (node != null) {
      // Bail out if we found a hit closer than the current node
      if (ray.maxDistance < tmin[0]) {
        break;
      }

      if (!node.isLeaf) {
        // Process kd-tree interior node
        Stats.KDTREE_INTERSECTIONP_TRAVERSED_INTERIOR_NODE(node);

        // Compute parametric distance along ray to split plane
        int axis = node.splitAxis;
        double tplane = (node.splitPos - ray.origin[axis]) * invDir[axis];

        // Get node children pointers for ray
        int firstChildIndex = 0;
        _KdAccelNode firstChild;
        int secondChildIndex = 0;
        _KdAccelNode secondChild;

        bool belowFirst = (ray.origin[axis] < node.splitPos) ||
                          (ray.origin[axis] == node.splitPos &&
                           ray.direction[axis] <= 0);

        if (belowFirst) {
          firstChildIndex = nodeIndex + 1;
          firstChild = nodes[firstChildIndex];
          secondChildIndex = node.aboveChild;
          secondChild = nodes[secondChildIndex];
        } else {
          firstChildIndex = node.aboveChild;
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
          // Enqueue secondChild in todo list
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
        Stats.KDTREE_INTERSECTIONP_TRAVERSED_LEAF_NODE(node, node.nPrimitives);

        int nPrimitives = node.nPrimitives;
        if (nPrimitives == 1) {
          Primitive prim = primitives[node.primitives];
          // Check one primitive inside leaf node
          Stats.KDTREE_INTERSECTIONP_PRIMITIVE_TEST(prim);
          if (prim.intersectP(ray)) {
            Stats.KDTREE_INTERSECTION_HIT(prim);
            return true;
          }
        } else if (nPrimitives > 1) {
          List<int> prims = node.primitives;
          for (int i = 0; i < nPrimitives; ++i) {
            Primitive prim = primitives[prims[i]];
            Stats.KDTREE_INTERSECTIONP_PRIMITIVE_TEST(prim);
            // Check one primitive inside leaf node
            if (prim.intersectP(ray)) {
              Stats.KDTREE_INTERSECTIONP_HIT(prim);
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

    Stats.KDTREE_INTERSECTIONP_MISSED();
    return false;
  }

  void _buildTree(int nodeNum, BBox nodeBounds, List<BBox> allPrimBounds,
                  Uint32List primNums, int depth,
                  List<List<_BoundEdge>> edges, Uint32List prims0,
                  Uint32List prims1, [int badRefines = 0]) {
    assert(nodeNum == nextFreeNode);

    final int nPrimitives = primNums.length;

    // Get next free node from _nodes_ array
    if (nodes == null || nextFreeNode == nodes.length) {
      final int numNodes = (nodes == null) ? 0 : nodes.length;

      int nAlloc = Math.max(2 * numNodes, 512);
      List<_KdAccelNode> n = new List<_KdAccelNode>(nAlloc);

      for (int i = 0; i < numNodes; ++i) {
        n[i] = nodes[i];
      }

      for (int i = numNodes; i < nAlloc; ++i) {
        n[i] = new _KdAccelNode();
      }

      nodes = n;
    }

    ++nextFreeNode;

    // Initialize leaf node if termination criteria met
    if (nPrimitives <= maxPrims || depth == 0) {
      Stats.KDTREE_CREATED_LEAF(nPrimitives, maxDepth - depth);
      nodes[nodeNum].initLeaf(primNums);
      return;
    }

    // Initialize interior node and continue recursion

    // Choose split axis position for interior node
    int bestAxis = -1;
    int bestOffset = -1;
    double bestCost = INFINITY;
    double oldCost = isectCost * nPrimitives.toDouble();
    Vector d = nodeBounds.pMax - nodeBounds.pMin;
    double totalSA = 2.0 * (d.x * d.y + d.x * d.z + d.y * d.z);
    double invTotalSA = 1.0 / totalSA;

    // Choose which axis to split along
    int axis = (d.x > d.y && d.x > d.z) ? 0 : (d.y > d.z) ? 1 : 2;

    const List<int> other_axis1 = const [1, 2, 0];
    const List<int> other_axis2 = const [2, 0, 1];

    int retries = 0;

    //retrySplit:
    while (true) {
      // Initialize edges for axis
      edges[axis] = new List<_BoundEdge>(nPrimitives * 2);
      for (int i = 0, j = 0; i < nPrimitives; ++i) {
        int pn = primNums[i];
        BBox bbox = allPrimBounds[pn];
        edges[axis][j++] = new _BoundEdge(bbox.pMin[axis], pn, true);
        edges[axis][j++] = new _BoundEdge(bbox.pMax[axis], pn, false);
      }

      edges[axis].sort((a, b) => a < b ? -1 : 1);

      // Compute cost of all splits for axis to find best
      int nBelow = 0;
      int nAbove = nPrimitives;
      for (int i = 0, len = 2 * nPrimitives; i < len; ++i) {
        if (edges[axis][i].type == _BoundEdge.END) {
          --nAbove;
        }

        double edget = edges[axis][i].t;
        if (edget > nodeBounds.pMin[axis] && edget < nodeBounds.pMax[axis]) {
          // Compute cost for split at i'th edge
          int otherAxis0 = other_axis1[axis];
          int otherAxis1 = other_axis2[axis];

          double belowSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                (edget - nodeBounds.pMin[axis]) *
                                (d[otherAxis0] + d[otherAxis1]));

          double aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                (nodeBounds.pMax[axis] - edget) *
                                (d[otherAxis0] + d[otherAxis1]));

          double pBelow = belowSA * invTotalSA;
          double pAbove = aboveSA * invTotalSA;
          double eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0.0;

          double cost = traversalCost + isectCost * (1.0 - eb) *
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
        axis = other_axis1[axis];
        continue;
      }

      break;
    }

    if (bestCost > oldCost) {
      ++badRefines;
    }

    if ((bestCost > 4.0 * oldCost && nPrimitives < 16) || bestAxis == -1 ||
        badRefines == 3) {
      Stats.KDTREE_CREATED_LEAF(nPrimitives, maxDepth - depth);
      nodes[nodeNum].initLeaf(primNums);
      return;
    }

    // Classify primitives with respect to split
    int n0 = 0;
    for (int i = 0; i < bestOffset; ++i) {
      if (edges[bestAxis][i].type == _BoundEdge.START) {
        prims0[n0++] = edges[bestAxis][i].primNum;
      }
    }

    int n1 = 0;
    for (int i = bestOffset + 1; i < 2 * nPrimitives; ++i) {
      if (edges[bestAxis][i].type == _BoundEdge.END) {
        prims1[n1++] = edges[bestAxis][i].primNum;
      }
    }

    // Recursively initialize children nodes
    double tsplit = edges[bestAxis][bestOffset].t;
    Stats.KDTREE_CREATED_INTERIOR_NODE(bestAxis, tsplit);

    BBox bounds0 = new BBox.from(nodeBounds);
    bounds0.pMax[bestAxis] = tsplit;

    BBox bounds1 = new BBox.from(nodeBounds);
    bounds1.pMin[bestAxis] = tsplit;

    // primNums = prims0, numPrimitives = n0
    Uint32List pn = new Uint32List.view(prims0.buffer, prims0.offsetInBytes,
                                        n0);

    // prims1 + nPrimitives
    Uint32List p1 = new Uint32List.view(prims1.buffer, prims1.offsetInBytes +
                                        (nPrimitives * 4));

    _buildTree(nodeNum + 1, bounds0, allPrimBounds, pn, depth - 1, edges,
               prims0, p1, badRefines);

    int aboveChild = nextFreeNode;
    nodes[nodeNum].initInterior(bestAxis, aboveChild, tsplit);

    pn = new Uint32List.view(prims1.buffer, prims1.offsetInBytes, n1);

    _buildTree(aboveChild, bounds1, allPrimBounds, pn, depth - 1, edges, prims0,
               p1, badRefines);
  }

  static KdTreeAccel Create(List<Primitive> prims, ParamSet ps) {
    int isectCost = ps.findOneInt('intersectcost', 80);
    int travCost = ps.findOneInt('traversalcost', 1);
    double emptyBonus = ps.findOneFloat('emptybonus', 0.5);
    int maxPrims = ps.findOneInt('maxprims', 1);
    int maxDepth = ps.findOneInt('maxdepth', -1);
    return new KdTreeAccel(prims, isectCost, travCost, emptyBonus, maxPrims,
                           maxDepth);
  }

  int isectCost;
  int traversalCost;
  int maxPrims;
  int maxDepth;
  double emptyBonus;
  List<Primitive> primitives = [];
  List<_KdAccelNode> nodes;
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
  void initLeaf(Uint32List primNums) {
    final int np = primNums.length;
    flags = 3;
    flags |= (np << 2);
    // Store primitive ids for leaf node
    if (np == 0) {
      primitives = 0;
    } else if (np == 1) {
      primitives = primNums[0];
    } else {
      primitives = new Uint32List.fromList(primNums);
    }
  }

  void initInterior(int axis, int aboveChild, double s) {
    split = s;
    flags = axis;
    flags |= (aboveChild << 2);
  }

  double get splitPos => split;

  int get nPrimitives => flags >> 2;

  int get aboveChild => flags >> 2;

  int get splitAxis => flags & 3;

  bool get isLeaf => (flags & 3) == 3;

  /// Bits 1-2: 0: x-axis leaf node; 1 y-axis leaf node; 2 z-axis leaf node;
  ///           3: interior node.
  /// Bits 3-n: the number of primitives stored by the node.
  int flags = 0;
  /// Split distance for interior nodes.
  double split;
  /// Stores an int if [nPrimitives] is 1, otherwise a [Uint32List].
  var primitives;
}


class _BoundEdge {
  _BoundEdge([this.t = 0.0, int primNum = 0, bool starting = false])
      : type_primNum = (starting ? START : END) | (primNum << 1);

  bool operator <(_BoundEdge e) {
    if (t == e.t) {
      return type > e.type;
    } else {
      return t < e.t;
    }
  }

  int get type => type_primNum & 1;

  int get primNum => type_primNum >> 1;

  static const int START = 1;
  static const int END = 0;

  double t;
  // type and primNum merged
  int type_primNum;
}
