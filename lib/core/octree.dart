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

class Octree {
  Octree(BBox b, [this.maxDepth = 16])
      : bound = new BBox.from(b),
        root = new _OctreeNode();

  void add(dataItem, BBox dataBound) {
    _add(root, bound, dataItem, dataBound,
               Vector.DistanceSquared(dataBound.pMin, dataBound.pMax));
  }

  void lookup(Point p, process) {
    if (!bound.inside(p)) {
      return;
    }

    _lookup(root, bound, p, process);
  }

  void _add(_OctreeNode node, BBox nodeBound, dataItem, BBox dataBound,
                  double diag2, [int depth = 0]) {
    // Possibly add data item to current octree node
    if (depth == maxDepth ||
        Vector.DistanceSquared(nodeBound.pMin, nodeBound.pMax) < diag2) {
      node.data.add(dataItem);
      return;
    }

    // Otherwise add data item to octree children
    Point pMid = nodeBound.pMin * 0.5 + nodeBound.pMax * 0.5;

    // Determine which children the item overlaps
    List<bool> x = [dataBound.pMin.x <= pMid.x, dataBound.pMax.x > pMid.x];
    List<bool> y = [dataBound.pMin.y <= pMid.y, dataBound.pMax.y > pMid.y];
    List<bool> z = [dataBound.pMin.z <= pMid.z, dataBound.pMax.z > pMid.z];
    List<bool> over = [(x[0] && y[0] && z[0]), (x[0] && y[0] && z[1]),
                       (x[0] && y[1] && z[0]), (x[0] && y[1] && z[1]),
                       (x[1] && y[0] && z[0]), (x[1] && y[0] && z[1]),
                       (x[1] && y[1] && z[0]), (x[1] && y[1] && z[1])];

    for (int child = 0; child < 8; ++child) {
      if (!over[child]) {
        continue;
      }

      // Allocate octree node if needed and continue recursive traversal
      if (node.children[child] == null) {
        node.children[child] = new _OctreeNode();
      }

      BBox childBound = OctreeChildBound(child, nodeBound, pMid);

      _add(node.children[child], childBound,
           dataItem, dataBound, diag2, depth + 1);
     }
  }

  bool _lookup(_OctreeNode node, BBox nodeBound, Point p, process) {
    for (int i = 0; i < node.data.length; ++i) {
      if (!process(node.data[i])) {
        return false;
      }
    }

    // Determine which octree child node _p_ is inside
    Point pMid = nodeBound.pMin * 0.5 + nodeBound.pMax * 0.5;

    int child = (p.x > pMid.x ? 4 : 0) + (p.y > pMid.y ? 2 : 0) +
                (p.z > pMid.z ? 1 : 0);

    if (node.children[child] == null) {
      return true;
    }

    BBox childBound = OctreeChildBound(child, nodeBound, pMid);

    return _lookup(node.children[child], childBound, p, process);
  }

  static BBox OctreeChildBound(int child, BBox nodeBound, Point pMid) {
     BBox childBound = new BBox();
     childBound.pMin.x = (child & 4) != 0 ? pMid.x : nodeBound.pMin.x;
     childBound.pMax.x = (child & 4) != 0 ? nodeBound.pMax.x : pMid.x;
     childBound.pMin.y = (child & 2) != 0 ? pMid.y : nodeBound.pMin.y;
     childBound.pMax.y = (child & 2) != 0 ? nodeBound.pMax.y : pMid.y;
     childBound.pMin.z = (child & 1) != 0 ? pMid.z : nodeBound.pMin.z;
     childBound.pMax.z = (child & 1) != 0 ? nodeBound.pMax.z : pMid.z;
     return childBound;
   }

  int maxDepth;
  BBox bound;
  _OctreeNode root;
}


class _OctreeNode {
  final List<_OctreeNode> children;
  final List data;

  _OctreeNode()
      : children = new List<_OctreeNode>(8),
        data = [];
}
