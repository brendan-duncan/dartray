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
 * A [KdTree] is a special case of Binary Space Partition, useful for
 * applications such as Nearest Neighbor searches.  You can use it, for
 * example, to find a set of points that are within some distance of a
 * point.
 */
class KdTree {
  KdTree(List data) {
    nNodes = data.length;
    nextFreeNode = 1;
    nodes = new List<_KdNode>(nNodes);
    nodeData = new List(nNodes);
    List<int> buildNodes = new List<int>(nNodes);
    for (int i = 0; i < nNodes; ++i) {
      buildNodes[i] = i;
    }
    // Begin the KdTree building process
    _recursiveBuild(0, 0, nNodes, data, buildNodes);
  }

  void lookup(Point p, process, List<double> maxDistSquared) {
    _lookup(0, p, process, maxDistSquared);
  }

  void _recursiveBuild(int nodeNum, int start, int end, List data,
                       List<int> buildNodes) {
    // Create leaf node of kd-tree if we've reached the bottom
    if (start + 1 == end) {
      nodes[nodeNum] = new _KdNode.leaf();
      nodeData[nodeNum] = data[buildNodes[start]];
      return;
    }

    // Choose split direction and partition data

    // Compute bounds of data from start to end
    BBox bound = new BBox();
    for (int i = start; i < end; ++i) {
      bound = BBox.UnionPoint(bound, data[buildNodes[i]].p);
    }

    int splitAxis = bound.maximumExtent();
    int splitPos = (start + end) ~/ 2;

    nth_element(buildNodes, start, splitPos, end,
                new _CompareNode(data, splitAxis));

    // Allocate kd-tree node and continue recursively
    nodes[nodeNum] = new _KdNode(data[buildNodes[splitPos]].p[splitAxis],
                                 splitAxis);
    nodeData[nodeNum] = data[buildNodes[splitPos]];

    if (start < splitPos) {
      nodes[nodeNum].hasLeftChild = true;
      int childNum = nextFreeNode++;
      _recursiveBuild(childNum, start, splitPos, data, buildNodes);
    }

    if (splitPos + 1 < end) {
      nodes[nodeNum].rightChild = nextFreeNode++;
      _recursiveBuild(nodes[nodeNum].rightChild, splitPos + 1,
                      end, data, buildNodes);
    }
  }

  void _lookup(int nodeNum, Point p, process, List<double> maxDistSquared) {
    _KdNode node = nodes[nodeNum];
    // Process kd-tree node's children
    int axis = node.splitAxis;
    if (axis != 3) {
      double dist2 = (p[axis] - node.splitPos) * (p[axis] - node.splitPos);
      if (p[axis] <= node.splitPos) {
        if (node.hasLeftChild) {
          _lookup(nodeNum + 1, p, process, maxDistSquared);
        }
        if (dist2 < maxDistSquared[0] && node.rightChild < nNodes) {
          _lookup(node.rightChild, p, process, maxDistSquared);
        }
      } else {
        if (node.rightChild < nNodes) {
          _lookup(node.rightChild, p, process, maxDistSquared);
        }
        if (dist2 < maxDistSquared[0] && node.hasLeftChild) {
          _lookup(nodeNum + 1, p, process, maxDistSquared);
        }
      }
    }

    // Hand kd-tree node to processing function
    double dist2 = Vector.DistanceSquared(nodeData[nodeNum].p, p);
    if (dist2 < maxDistSquared[0]) {
      process(p, nodeData[nodeNum], dist2, maxDistSquared);
    }
  }

  List<_KdNode> nodes;
  List nodeData;
  int nNodes;
  int nextFreeNode;
}

class _CompareNode {
  _CompareNode(this.data, this.axis);

  bool call(d1, d2) {
    return (data[d1].p[axis] == data[d2].p[axis]) ?
           (data[d1].hashCode < data[d2].hashCode) :
           (data[d1].p[axis] < data[d2].p[axis]);
  }

  List data;
  int axis;
}

class _KdNode {
  _KdNode(this.splitPos, this.splitAxis)
      : rightChild = (1 << 29) - 1,
        hasLeftChild = false;

  _KdNode.leaf()
      : splitAxis = 3,
        rightChild = (1 << 29) - 1,
        hasLeftChild = false;

  double splitPos;
  int splitAxis;
  bool hasLeftChild;
  int rightChild;
}
