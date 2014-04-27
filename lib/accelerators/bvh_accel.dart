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
 * Bounding volume hierarchy ray accelerator.
 */
class BVHAccel extends Aggregate {
  static const int SPLIT_MIDDLE = 0;
  static const int SPLIT_EQUAL_COUNTS = 1;
  static const int SPLIT_SAH = 2;

  BVHAccel(List<Primitive> p, [int maxPrims = 1,
           this.splitMethod = SPLIT_SAH]) {
    LogInfo('Building BVH Acceleration Structures.');
    maxPrimsInNode = Math.min(255, maxPrims);
    for (int i = 0; i < p.length; ++i) {
      p[i].fullyRefine(primitives);
    }
    LogInfo('BVH: ${primitives.length} Primitives');

    if (primitives.isEmpty) {
      nodes = null;
      return;
    }

    // Build BVH from primitives
    Stats.BVH_STARTED_CONSTRUCTION(this, primitives.length);

    // Initialize _buildData_ array for primitives
    List<_BVHPrimitiveInfo> buildData =
        new List<_BVHPrimitiveInfo>(primitives.length);

    for (int i = 0; i < primitives.length; ++i) {
      BBox bbox = primitives[i].worldBound();
      buildData[i] = new _BVHPrimitiveInfo(i, bbox);
    }

    // Recursively build BVH tree for primitives
    List<int> totalNodes = [0];
    List<Primitive> orderedPrims = [];

    _BVHBuildNode root = _recursiveBuild(buildData, 0,
                                         primitives.length,
                                         totalNodes,
                                         orderedPrims);

    primitives = orderedPrims;

    LogInfo('BVH created with ${totalNodes[0]} nodes for '
            '${primitives.length} primitives');

    // Compute representation of depth-first traversal of BVH tree
    nodes = new List<_LinearBVHNode>(totalNodes[0]);
    for (int i = 0, len = totalNodes[0]; i < len; ++i) {
      nodes[i] = new _LinearBVHNode();
    }

    List<int> offset = [0];
    _flattenBVHTree(root, offset);
    assert(offset[0] == totalNodes[0]);
    Stats.BVH_FINISHED_CONSTRUCTION(this);
  }

  BBox worldBound() {
    return nodes != null ? nodes[0].bounds : new BBox();
  }

  bool canIntersect() {
    return true;
  }

  bool intersect(Ray ray, Intersection isect) {
    if (nodes == null) {
      return false;
    }

    Stats.BVH_INTERSECTION_STARTED(this, ray);

    bool hit = false;
    final Vector invDir = new Vector(1.0 / ray.direction.x,
                                     1.0 / ray.direction.y,
                                     1.0 / ray.direction.z);

    final List<int> dirIsNeg = [invDir.x < 0 ? 1 : 0,
                                invDir.y < 0 ? 1 : 0,
                                invDir.z < 0 ? 1 : 0];

    // Follow ray through BVH nodes to find primitive intersections
    int todoOffset = 0;
    int nodeNum = 0;
    final Uint32List todo = new Uint32List(64);

    while (true) {
      _LinearBVHNode node = nodes[nodeNum];
      // Check ray against BVH node
      if (_intersectP(node.bounds, ray, invDir, dirIsNeg)) {
        if (node.nPrimitives > 0) {
          // Intersect ray with primitives in leaf BVH node
          Stats.BVH_INTERSECTION_TRAVERSED_LEAF_NODE(node);
          for (int i = 0; i < node.nPrimitives; ++i) {
            Stats.BVH_INTERSECTION_PRIMITIVE_TEST(primitives[node.offset + i]);
            if (primitives[node.offset + i].intersect(ray, isect)) {
              Stats.BVH_INTERSECTION_PRIMITIVE_HIT(primitives[node.offset + i]);
              hit = true;
            } else {
              Stats.BVH_INTERSECTION_PRIMITIVE_MISSED(primitives[node.offset + i]);
            }
          }

          if (todoOffset == 0) {
            break;
          }

          nodeNum = todo[--todoOffset];
        } else {
          // Put far BVH node on _todo_ stack, advance to near node
          Stats.BVH_INTERSECTION_TRAVERSED_INTERIOR_NODE(node);
          if (dirIsNeg[node.axis] != 0) {
            todo[todoOffset++] = nodeNum + 1;
            nodeNum = node.offset;
          } else {
            todo[todoOffset++] = node.offset;
            nodeNum = nodeNum + 1;
          }
        }
      } else {
        if (todoOffset == 0) {
          break;
        }
        nodeNum = todo[--todoOffset];
      }
    }

    Stats.BVH_INTERSECTION_FINISHED();
    return hit;
  }

  bool intersectP(Ray ray) {
    if (nodes == null) {
      return false;
    }
    Stats.BVH_INTERSECTIONP_STARTED(this, ray);

    final Vector invDir = new Vector(1.0 / ray.direction.x,
                                     1.0 / ray.direction.y,
                                     1.0 / ray.direction.z);

    final List<int> dirIsNeg = [invDir.x < 0 ? 1 : 0,
                                invDir.y < 0 ? 1 : 0,
                                invDir.z < 0 ? 1 : 0];

    final Uint32List todo = new Uint32List(64);
    int todoOffset = 0;
    int nodeNum = 0;

    while (true) {
      _LinearBVHNode node = nodes[nodeNum];
      if (_intersectP(node.bounds, ray, invDir, dirIsNeg)) {
        // Process BVH node _node_ for traversal
        if (node.nPrimitives > 0) {
          Stats.BVH_INTERSECTIONP_TRAVERSED_LEAF_NODE(node);
          for (int i = 0; i < node.nPrimitives; ++i) {
            Stats.BVH_INTERSECTIONP_PRIMITIVE_TEST(primitives[node.offset + i]);
            if (primitives[node.offset + i].intersectP(ray)) {
              Stats.BVH_INTERSECTIONP_PRIMITIVE_HIT(primitives[node.offset + i]);
              return true;
            } else {
              Stats.BVH_INTERSECTIONP_PRIMITIVE_MISSED(primitives[node.offset + i]);
            }
          }
          if (todoOffset == 0) {
            break;
          }

          nodeNum = todo[--todoOffset];
        } else {
          Stats.BVH_INTERSECTIONP_TRAVERSED_INTERIOR_NODE(node);
          if (dirIsNeg[node.axis] != 0) {
            // second child first
            todo[todoOffset++] = nodeNum + 1;
            nodeNum = node.offset;
          } else {
            todo[todoOffset++] = node.offset;
            nodeNum = nodeNum + 1;
          }
        }
      } else {
        if (todoOffset == 0) {
          break;
        }
        nodeNum = todo[--todoOffset];
      }
    }

    Stats.BVH_INTERSECTIONP_FINISHED();
    return false;
  }

  _BVHBuildNode _recursiveBuild(List<_BVHPrimitiveInfo> buildData,
                                int start,
                                int end,
                                List<int> totalNodes,
                                List<Primitive> orderedPrims) {
    assert(start != end);
    totalNodes[0]++;

    _BVHBuildNode node = new _BVHBuildNode();
    // Compute bounds of all primitives in BVH node
    BBox bbox = new BBox();
    for (int i = start; i < end; ++i) {
      bbox = BBox.Union(bbox, buildData[i].bounds);
    }

    int nPrimitives = end - start;

    if (nPrimitives == 1) {
      // Create leaf _BVHBuildNode_
      int firstPrimOffset = orderedPrims.length;
      for (int i = start; i < end; ++i) {
        int primNum = buildData[i].primitiveNumber;
        orderedPrims.add(primitives[primNum]);
      }
      node.initLeaf(firstPrimOffset, nPrimitives, bbox);
    } else {
      // Compute bound of primitive centroids, choose split dimension _dim_
      BBox centroidBounds = new BBox();
      for (int i = start; i < end; ++i) {
        centroidBounds = BBox.UnionPoint(centroidBounds, buildData[i].centroid);
      }

      int dim = centroidBounds.maximumExtent();

      // Partition primitives into two sets and build children
      int mid = (start + end) ~/ 2;

      if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
        // Create leaf _BVHBuildNode_
        int firstPrimOffset = orderedPrims.length;
        for (int i = start; i < end; ++i) {
          int primNum = buildData[i].primitiveNumber;
          orderedPrims.add(primitives[primNum]);
        }
        node.initLeaf(firstPrimOffset, nPrimitives, bbox);
        return node;
      }

      bool ComparePoints(_BVHPrimitiveInfo a, _BVHPrimitiveInfo b) {
        return a.centroid[dim] < b.centroid[dim];
      }

      // Partition primitives based on _splitMethod_
      switch (splitMethod) {
        case SPLIT_MIDDLE:
          // Partition primitives through node's midpoint
          double pmid = 0.5 * (centroidBounds.pMin[dim] +
                               centroidBounds.pMax[dim]);

          bool CompareToMid(_BVHPrimitiveInfo a) {
            return a.centroid[dim] < pmid;
          }

          mid = partition(buildData, CompareToMid, start, end);

          if (mid != start && mid != end) {
            // for lots of prims with large overlapping bounding boxes, this
            // may fail to partition; in that case don't break and fall through
            // to SPLIT_EQUAL_COUNTS
            break;
          }

          // Partition primitives into equally-sized subsets
          mid = (start + end) ~/ 2;

          nth_element(buildData, start, mid, end, ComparePoints);
          break;
        case SPLIT_EQUAL_COUNTS:
          // Partition primitives into equally-sized subsets
          mid = (start + end) ~/ 2;
          nth_element(buildData, start, mid, end, ComparePoints);
          break;
        case SPLIT_SAH:
        default:
          // Partition primitives using approximate SAH
          if (nPrimitives <= 4) {
            // Partition primitives into equally-sized subsets
            mid = (start + end) ~/ 2;
            nth_element(buildData, start, mid, end, ComparePoints);
          } else {
            // Allocate _BucketInfo_ for SAH partition buckets
            const int nBuckets = 12;

            List<_BVHBucketInfo> buckets = new List<_BVHBucketInfo>(nBuckets);
            for (int i = 0; i < nBuckets; ++i) {
              buckets[i] = new _BVHBucketInfo();
            }

            // Initialize _BucketInfo_ for SAH partition buckets
            for (int i = start; i < end; ++i) {
              int b = (nBuckets *
                ((buildData[i].centroid[dim] - centroidBounds.pMin[dim]) /
                (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]))).toInt();
              if (b == nBuckets) {
                b = nBuckets - 1;
              }

              assert(b >= 0 && b < nBuckets);
              buckets[b].count++;
              buckets[b].bounds = BBox.Union(buckets[b].bounds,
                                             buildData[i].bounds);
            }

            BBox b0 = new BBox();
            BBox b1 = new BBox();

            // Compute costs for splitting after each bucket
            List<double> cost = new Float32List(nBuckets - 1);

            for (int i = 0; i < nBuckets - 1; ++i) {
              b0.reset();
              b1.reset();
              int count0 = 0, count1 = 0;

              for (int j = 0; j <= i; ++j) {
                b0 = BBox.Union(b0, buckets[j].bounds);
                count0 += buckets[j].count;
              }

              for (int j = i + 1; j < nBuckets; ++j) {
                b1 = BBox.Union(b1, buckets[j].bounds);
                count1 += buckets[j].count;
              }

              cost[i] = 0.125 + (count0 * b0.surfaceArea() +
                                count1 * b1.surfaceArea()) / bbox.surfaceArea();
            }

            // Find bucket to split at that minimizes SAH metric
            double minCost = cost[0];
            int minCostSplit = 0;
            for (int i = 1; i < nBuckets - 1; ++i) {
              if (cost[i] < minCost) {
                minCost = cost[i];
                minCostSplit = i;
              }
            }

            // Either create leaf or split primitives at selected SAH bucket
            if (nPrimitives > maxPrimsInNode || minCost < nPrimitives) {
              bool CompareToBucket(_BVHPrimitiveInfo p) {
                int b = (nBuckets *
                    ((p.centroid[dim] - centroidBounds.pMin[dim]) /
                     (centroidBounds.pMax[dim] -
                      centroidBounds.pMin[dim]))).floor();

                if (b == nBuckets) {
                  b = nBuckets - 1;
                }

                assert(b >= 0 && b < nBuckets);
                return b <= minCostSplit;
              }

              mid = partition(buildData, CompareToBucket, start, end);
            } else {
              // Create leaf _BVHBuildNode_
              int firstPrimOffset = orderedPrims.length;
              for (int i = start; i < end; ++i) {
                int primNum = buildData[i].primitiveNumber;
                orderedPrims.add(primitives[primNum]);
              }
              node.initLeaf(firstPrimOffset, nPrimitives, bbox);
              return node;
            }
          }
          break;
      }

      _BVHBuildNode c2 = _recursiveBuild(buildData, mid, end,
                                         totalNodes, orderedPrims);

      _BVHBuildNode c1 = _recursiveBuild(buildData, start, mid,
                                         totalNodes, orderedPrims);

      node.initInterior(dim, c1, c2);
    }

    return node;
  }

  int _flattenBVHTree(_BVHBuildNode node, List<int> offset) {
    _LinearBVHNode linearNode = nodes[offset[0]];
    linearNode.bounds = node.bounds;//new BBox.from(node.bounds);
    int myOffset = offset[0]++;

    if (node.nPrimitives > 0) {
      assert(node.children[0] == null && node.children[1] == null);
      linearNode.offset = node.firstPrimOffset;
      linearNode.nPrimitives = node.nPrimitives;
    } else {
      // Creater interior flattened BVH node
      linearNode.axis = node.splitAxis;
      linearNode.nPrimitives = 0;
      _flattenBVHTree(node.children[0], offset);
      linearNode.offset = _flattenBVHTree(node.children[1], offset);
    }

    return myOffset;
  }

  static _intersectP(BBox bounds, Ray ray, Vector invDir, List<int> dirIsNeg) {
    // Check for ray intersection against $x$ and $y$ slabs
    double tmin = (bounds[dirIsNeg[0]].x - ray.origin.x) * invDir.x;
    double tmax = (bounds[1 - dirIsNeg[0]].x - ray.origin.x) * invDir.x;
    double tymin = (bounds[dirIsNeg[1]].y - ray.origin.y) * invDir.y;
    double tymax = (bounds[1 - dirIsNeg[1]].y - ray.origin.y) * invDir.y;

    if ((tmin > tymax) || (tymin > tmax)) {
      return false;
    }

    if (tymin > tmin) {
      tmin = tymin;
    }
    if (tymax < tmax) {
      tmax = tymax;
    }

    // Check for ray intersection against z slab
    double tzmin = (bounds[dirIsNeg[2]].z - ray.origin.z) * invDir.z;
    double tzmax = (bounds[1 - dirIsNeg[2]].z - ray.origin.z) * invDir.z;
    if ((tmin > tzmax) || (tzmin > tmax)) {
      return false;
    }

    if (tzmin > tmin) {
      tmin = tzmin;
    }
    if (tzmax < tmax) {
      tmax = tzmax;
    }

    return (tmin < ray.maxDistance) && (tmax > ray.minDistance);
  }

  static BVHAccel Create(List<Primitive> prims, ParamSet ps) {
    String splitMethod = ps.findOneString("splitmethod", "sah");
    int maxPrimsInNode = ps.findOneInt("maxnodeprims", 4);
    int sm = (splitMethod == 'sah') ? SPLIT_SAH :
             (splitMethod == 'middle') ? SPLIT_MIDDLE :
             (splitMethod == 'equal') ? SPLIT_EQUAL_COUNTS :
             SPLIT_SAH;
    return new BVHAccel(prims, maxPrimsInNode, sm);
  }

  int maxPrimsInNode;
  int splitMethod;
  List<Primitive> primitives = [];
  List<_LinearBVHNode> nodes;
}

class _BVHPrimitiveInfo {
  _BVHPrimitiveInfo([this.primitiveNumber = 0, this.bounds]) {
    if (bounds == null) {
      bounds = new BBox();
    }
    centroid = bounds.center;
  }

  int primitiveNumber;
  Point centroid;
  BBox bounds;
}

class _BVHBucketInfo {
  int count = 0;
  BBox bounds = new BBox();
}

class _BVHBuildNode {
  _BVHBuildNode() :
    bounds = new BBox();

  void initLeaf(int first, int n, BBox b) {
    firstPrimOffset = first;
    nPrimitives = n;
    bounds.copy(b);
  }

  void initInterior(int axis, _BVHBuildNode c0, _BVHBuildNode c1) {
    children[0] = c0;
    children[1] = c1;
    bounds = BBox.Union(c0.bounds, c1.bounds);
    splitAxis = axis;
    nPrimitives = 0;
  }

  BBox bounds;
  List<_BVHBuildNode> children = [null, null];
  int splitAxis;
  int firstPrimOffset;
  int nPrimitives;
}

class _LinearBVHNode {
  BBox bounds;
  int offset; // primitivesOffset / secondChildOffset
  int nPrimitives; // 0 . interior node
  int axis;
}
