/****************************************************************************
 *  Copyright (C) 2014 by Brendan Duncan.                                   *
 *                                                                          *
 *  This file is part of DartRay.                                           *
 *                                                                          *
 *  Licensed under the Apache License, Version 2.0 (the "License");         *
 *  you may not use this file except in compliance with the License.        *
 *  You may obtain a copy of the License at                                 *
 *                                                                          *
 *  http://www.apache.org/licenses/LICENSE-2.0                              *
 *                                                                          *
 *  Unless required by applicable law or agreed to in writing, software     *
 *  distributed under the License is distributed on an "AS IS" BASIS,       *
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 *  See the License for the specific language governing permissions and     *
 *  limitations under the License.                                          *
 *                                                                          *
 *   This project is based on PBRT v2 ; see http://www.pbrt.org             *
 *   pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.*
 ****************************************************************************/
part of core;


/**
 * A [KdTree] is a special case of Binary Space Partition, useful for
 * applications such as Nearest Neighbor searches.  You can use it, for
 * example, to find a set of points that are within some distance of a
 * point.
 */
class KdTree {
  static const int X_AXIS = 0;
  static const int Y_AXIS = 1;
  static const int Z_AXIS = 2;

  KdTree(List data) {
    if (data != null) {
      setup(data);
    }
  }

  /**
   * Reset the tree.
   */
  void clear() {
    tree = null;
    data = null;
  }

  /**
   * Does the tree have data?
   */
  bool get isValid => data != null;

  /**
   * Constructs the kd-tree from the given [data].
   */
  void setup(List data) {
    this.data = data;
    if (data == null || data.length == 0) {
      return;
    }

    tree = new List<_KdTreeNode>(data.length);
    tree[0] = new _KdTreeNode();

    final int num = data.length;
    Point p = data[0].p;
    bounds = new BBox(p);
    for (int i = 1; i < num; ++i) {
      tree[i] = new _KdTreeNode();
      bounds.unionPoint(data[i].p);
      tree[i].index = i;
    }

    _sortAndSubdivide(0, data.length - 1, X_AXIS + Y_AXIS + Z_AXIS);
  }

  void lookup(Vector p, process, List<double> radius) {
    _lookup(p, 0, numPoints - 1, radius, process);
  }

  /**
   * Find the points in the tree that are within [radius] distance of [p].
   */
  int findNearest(Vector p, [double radius = 1.0e10, List<int> result,
                  List<double> distances, bool sorted = false]) {
    if (result == null) {
      result = new List<int>();
      sorted = false;
    }

    if (distances == null) {
      distances = new List<double>();
      sorted = false;
    }

    List nearest = [-1, INFINITY];
    _gather(p, 0, numPoints - 1, radius, result, distances, nearest);

    if (sorted) {
      _sort(result, distances);
    }

    return result.length;
  }

  /**
   * How many points in the kd-tree.
   */
  int get numPoints => data.length;


  Point _position(int treeIndex) {
    return data[tree[treeIndex].index].p;
  }

  int _axis(int treeIndex) {
    return tree[treeIndex].axis;
  }

  int _dataIndex(int treeIndex) {
    return tree[treeIndex].index;
  }

  void _sortAndSubdivide(int start, int end, int sorted) {
    if (end == start) {
      tree[start].axis = -1;
      return;
    } else if (end < start) {
      return;
    }

    // Go over the points, building the bounding box for them to find
    // the dominant axis that the points span.
    Point p = _position(start);
    BBox bounds = new BBox(p);
    for (int i = start + 1; i <= end; ++i) {
      p = _position(i);
      bounds.unionPoint(p);
    }

    int dominantAxis = bounds.maximumExtent();

    if (dominantAxis != sorted) {
      int len = end - start;
      if (len > 1000) {
        _quickSortRec(start, end, dominantAxis);
      }
      _insertSort(start, end, dominantAxis);
    }

    // Find the midpoint and set dominantAxis as its axis.
    int mid = (end + start) >> 1;

    tree[mid].axis = dominantAxis;

    // Now recurse to continue building the kd-tree.
    _sortAndSubdivide(start, mid - 1, dominantAxis);
    _sortAndSubdivide(mid + 1, end, dominantAxis);
  }

  void _quickSortRec(int left, int right, int d) {
    if (left < right) {
      _swap(((left + right) >> 1), left + 1);

      if (_position(left + 1)[d] > _position(right)[d]) {
        _swap(left + 1, right);
      }

      if (_position(left)[d] > _position(right)[d]) {
        _swap(left, right);
      }

      if (_position(left + 1)[d] > _position(left)[d]) {
        _swap(left + 1, left);
      }

      int j = left + 1;
      int k = right;
      while (j <= k) {
        for (j++; ((j <= right) && (_position(j)[d] < _position(left)[d])); j++);
        for (k--; ((k >= left) && (_position(k)[d] > _position(left)[d])); k--);
        if (j < k) {
          _swap(j, k);
        }
      }

      _swap(left, k);

      if (k - left > 10) {
        _quickSortRec(left, k - 1, d);
      }

      if (right - k > 10) {
        _quickSortRec(k + 1, right, d);
      }

      // leave the rest for insertSort
    }
  }

  void _insertSort(int start, int end, int d) {
    Vector tmpPos = new Vector();
    for (int k = end - 1; k >= start; k--) {
      int j = k + 1;
      _KdTreeNode tmp = tree[k];
      tmpPos = _position(k);
      while (tmpPos[d] > _position(j)[d]) {
        tree[j - 1] = tree[j];
        j++;
        if (j > end) {
          break;
        }
      }
      tree[j - 1] = tmp;
    }
  }

  void _swap(int a, int b) {
    _KdTreeNode tmp = tree[a];
    tree[a] = tree[b];
    tree[b] = tmp;
  }

  /**
   * Sort the results of findNearest.
   */
  void _sort(List<int> result, List<double> distances) {
    int len = result.length;
    for (int k = len - 2; k >= 0; k--) {
      int j = k + 1;
      double tmpD = distances[k];
      int tmpP = result[k];
      while (tmpD > distances[j]) {
        distances[j - 1]  = distances[j];
        result[j - 1] = result[j];
        j++;
        if (j >= len) {
          break;
        }
      }
      distances[j - 1]  = tmpD;
      result[j - 1] = tmpP;
    }
  }

  void _lookup(Point p, int start, int end, List<double> radius, process) {
    // find midpoint
    double radius2 = radius[0] * radius[0];
    int mid = (end + start) >> 1;
    int axis = tree[mid].axis;
    Point P = _position(mid);

    // check this element

    // find distance from pt
    Point delta = p - P;
    double dist2 = delta.lengthSquared();

    if (dist2 < radius2) {
      process(p, data[mid], dist2, radius);
    }

    // now go left & right if appropriate - if going left or
    // right goes out the current range, then don't go that way.

    if (delta[axis] < 0.0) {
      // on left - go left first
      if (P[axis] - radius[0] < P[axis]) {
        if (mid - 1 >= start) {
          _lookup(p, start, mid - 1, radius, process);
        }
      }

      if (p[axis] + radius[0] > P[axis]) {
        if (end >= mid + 1) {
          _lookup(p, mid + 1, end, radius, process);
        }
      }
    } else {
      // on right - go right first
      if (p[axis] + radius[0] > P[axis]) {
        if (end >= mid + 1) {
          _lookup(p, mid + 1, end, radius, process);
        }
      }

      if (p[axis] - radius[0] < P[axis]) {
        if (mid - 1 >= start) {
          _lookup(p, start, mid - 1, radius, process);
        }
      }
    }
  }

  void _gather(Point p, int start, int end, double radius,
              List<int> result, List<double> distances,
              List nearest) {
    // find midpoint
    double radius2 = radius * radius;
    int mid = (end + start) >> 1;
    int axis = tree[mid].axis;
    Point P = _position(mid);

    // check this element

    // find distance from pt
    Point delta = p - P;
    double dist2 = delta.lengthSquared();

    if (dist2 < radius2) {
      if (dist2 < nearest[1]) {
        nearest[0] = mid;
        nearest[1] = dist2;
      }
      result.add(mid);
      distances.add(dist2);
    }

    // now go left & right if appropriate - if going left or
    // right goes out the current range, then don't go that way.

    if (delta[axis] < 0.0) {
      // on left - go left first
      if (P[axis] - radius < P[axis]) {
        if (mid - 1 >= start) {
          _gather(p, start, mid - 1, radius, result, distances, nearest);
        }
      }

      if (p[axis] + radius > P[axis]) {
        if (end >= mid + 1) {
          _gather(p, mid + 1, end, radius, result, distances, nearest);
        }
      }
    } else {
      // on right - go right first
      if (p[axis] + radius > P[axis]) {
        if (end >= mid + 1) {
          _gather(p, mid + 1, end, radius, result, distances, nearest);
        }
      }

      if (p[axis] - radius < P[axis]) {
        if (mid - 1 >= start) {
          _gather(p, start, mid - 1, radius, result, distances, nearest);
        }
      }
    }
  }

  List data;
  BBox bounds;
  List<_KdTreeNode> tree;
}


/**
 * A node in a [KdTree].
 */
class _KdTreeNode {
  int index;
  int axis;

  _KdTreeNode([this.index = 0, this.axis = -1]);
}


/*class KdNode {
  void init(double p, int a) {
    splitPos = p;
    splitAxis = a;
    rightChild = (1 << 29) - 1;
    hasLeftChild = false;
  }

  void initLeaf() {
    splitAxis = 3;
    rightChild = (1 << 29) - 1;
    hasLeftChild = false;
  }

  double splitPos;
  int splitAxis;
  bool hasLeftChild;
  int rightChild;
}

class KdTree {
  KdTree(List data) {
    nNodes = data.length;
    nextFreeNode = 1;
    nodes = new List<KdNode>(nNodes);
    nodeData = new List(nNodes);
    List<int> buildNodes = new List<int>(nNodes);
    for (int i = 0; i < nNodes; ++i) {
      buildNodes[i] = i;
    }
    // Begin the KdTree building process
    _recursiveBuild(0, 0, nNodes, data, buildNodes, 0);
  }

  void lookup(Point p, process, List<double> maxDistSquared) {

  }

  static compareNode(int axis, d1, d2) =>
    d1.p[axis] == d2.p[axis] ? (d1 < d2) : d1.p[axis] < d2.p[axis];

  void _recursiveBuild(int nodeNum, int start, int end, List data,
                       List<int> buildNodes, int index) {
    // Create leaf node of kd-tree if we've reached the bottom
    if (start + 1 == end) {
      nodes[nodeNum].initLeaf();
      nodeData[nodeNum] = data[buildNodes[index + start]];
      return;
    }

    // Choose split direction and partition data

    // Compute bounds of data from _start_ to _end_
    BBox bound = new BBox();
    for (int i = start; i < end; ++i) {
      bound = BBox.Union(bound, data[index + buildNodes[i]].p);
    }

    int splitAxis = bound.maximumExtent();
    int splitPos = (start + end) ~/ 2;
    for (int li = start, gi = splitPos + 1; li < splitPos; ++li, ++gi) {

    }

    /*std::nth_element(&buildNodes[start], &buildNodes[splitPos],
                  &buildNodes[end], CompareNode(splitAxis));*/

    // Allocate kd-tree node and continue recursively
    nodes[nodeNum].init(data[index + buildNodes[splitPos]].p[splitAxis], splitAxis);
    nodeData[nodeNum] = data[index + buildNodes[splitPos]];

    if (start < splitPos) {
      nodes[nodeNum].hasLeftChild = true;
      int childNum = nextFreeNode++;
      _recursiveBuild(childNum, start, splitPos, data, buildNodes, index);
    }

    if (splitPos+1 < end) {
      nodes[nodeNum].rightChild = nextFreeNode++;
      _recursiveBuild(nodes[nodeNum].rightChild, splitPos + 1,
                      end, data, buildNodes, index);
    }
  }

  void _lookup(int nodeNum, Point p, process, List<double> maxDistSquared) {

  }

  List<KdNode> nodes;
  List nodeData;
  int nNodes;
  int nextFreeNode;
}*/
