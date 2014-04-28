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
 * Hierarchical uniform grid ray accelerator.
 */
class GridAccel extends Aggregate {
  GridAccel(List<Primitive> p, bool refineImmediately) {
    LogInfo('Building Hierarchical Grid Acceleration Structures.');
    Stats.GRID_STARTED_CONSTRUCTION(this, p.length);
    // Initialize _primitives_ with primitives for grid
    if (refineImmediately) {
      for (int i = 0; i < p.length; ++i) {
        p[i].fullyRefine(primitives);
      }
    } else {
      primitives = p;
    }

    // Compute bounds and choose grid resolution
    for (int i = 0; i < primitives.length; ++i) {
      bounds = BBox.Union(bounds, primitives[i].worldBound());
    }

    Vector delta = bounds.pMax - bounds.pMin;

    // Find _voxelsPerUnitDist_ for grid
    int maxAxis = bounds.maximumExtent();
    double invMaxWidth = 1.0 / delta[maxAxis];
    assert(invMaxWidth > 0.0);

    double cubeRoot = 3.0 * Math.pow(primitives.length, 1.0 / 3.0);
    double voxelsPerUnitDist = cubeRoot * invMaxWidth;
    for (int axis = 0; axis < 3; ++axis) {
      nVoxels[axis] = (delta[axis] * voxelsPerUnitDist).round();
      nVoxels[axis] = nVoxels[axis].clamp(1, 64);
    }
    Stats.GRID_BOUNDS_AND_RESOLUTION(bounds, nVoxels);

    // Compute voxel widths and allocate voxels
    for (int axis = 0; axis < 3; ++axis) {
      width[axis] = delta[axis] / nVoxels[axis];
      invWidth[axis] = (width[axis] == 0.0) ? 0.0 : 1.0 / width[axis];
    }

    int nv = nVoxels[0] * nVoxels[1] * nVoxels[2];
    voxels = new List<_Voxel>(nv);

    List<int> vmin = [0, 0, 0];
    List<int> vmax = [0, 0, 0];

    // Add primitives to grid voxels
    for (int i = 0; i < primitives.length; ++i) {
      Primitive prim = primitives[i];

      // Find voxel extent of primitive
      BBox pb = prim.worldBound();

      for (int axis = 0; axis < 3; ++axis) {
        vmin[axis] = posToVoxel(pb.pMin, axis);
        vmax[axis] = posToVoxel(pb.pMax, axis);
      }

      // Add primitive to overlapping voxels
      Stats.GRID_VOXELIZED_PRIMITIVE(vmin, vmax);
      for (int z = vmin[2]; z <= vmax[2]; ++z) {
        for (int y = vmin[1]; y <= vmax[1]; ++y) {
          for (int x = vmin[0]; x <= vmax[0]; ++x) {
            int o = offset(x, y, z);
            if (voxels[o] == null) {
              voxels[o] = new _Voxel(prim);
            } else {
              voxels[o].addPrimitive(prim);
            }
          }
        }
      }
    }
    Stats.GRID_FINISHED_CONSTRUCTION(this);
  }

  BBox worldBound() => new BBox.from(bounds);

  bool canIntersect() {
    return true;
  }

  bool intersect(Ray ray, Intersection isect) {
    Stats.GRID_INTERSECTION_TEST(this, ray);
    // Check ray against overall grid bounds
    List<double> rayT = [0.0];
    if (bounds.inside(ray.pointAt(ray.minDistance))) {
      rayT[0] = ray.minDistance;
    } else if (!bounds.intersectP(ray, rayT)) {
      Stats.GRID_RAY_MISSED_BOUNDS();
      return false;
    }

    Point gridIntersect = ray.pointAt(rayT[0]);

    // Set up 3D DDA for ray
    List<double> nextCrossingT = [0.0, 0.0, 0.0];
    List<double> deltaT = [0.0, 0.0, 0.0];
    List<int> step = [0, 0, 0];
    List<int> out = [0, 0, 0];
    List<int> pos = [0, 0, 0];

    for (int axis = 0; axis < 3; ++axis) {
      // Compute current voxel for axis
      pos[axis] = posToVoxel(gridIntersect, axis);

      if (ray.direction[axis] >= 0.0) {
        // Handle ray with positive direction for voxel stepping
        nextCrossingT[axis] = rayT[0] +
                 (voxelToPos(pos[axis] + 1, axis) - gridIntersect[axis]) /
                 ray.direction[axis];
        deltaT[axis] = width[axis] / ray.direction[axis];
        step[axis] = 1;
        out[axis] = nVoxels[axis];
      } else {
        // Handle ray with negative direction for voxel stepping
        nextCrossingT[axis] = rayT[0] +
                 (voxelToPos(pos[axis], axis) - gridIntersect[axis]) /
                 ray.direction[axis];
        deltaT[axis] = -width[axis] / ray.direction[axis];
        step[axis] = -1;
        out[axis] = -1;
      }
    }

    // Walk ray through voxel grid
    bool hitSomething = false;
    while (true) {
      // Check for intersection in current voxel and advance to next
      _Voxel voxel = voxels[offset(pos[0], pos[1], pos[2])];
      Stats.GRID_RAY_TRAVERSED_VOXEL(pos, voxel != null ? voxel.size() : 0);
      if (voxel != null) {
        if (voxel.intersect(ray, isect)) {
          hitSomething = true;
        }
      }

      // Advance to next voxel

      // Find _stepAxis_ for stepping to next voxel
      int bits = ((nextCrossingT[0] < nextCrossingT[1]) ? 4 : 0) +
                 ((nextCrossingT[0] < nextCrossingT[2]) ? 2 : 0) +
                 ((nextCrossingT[1] < nextCrossingT[2]) ? 1 : 0);
      const List<int> cmpToAxis = const [ 2, 1, 2, 1, 2, 2, 0, 0 ];
      int stepAxis = cmpToAxis[bits];

      if (ray.maxDistance < nextCrossingT[stepAxis]) {
        break;
      }

      pos[stepAxis] += step[stepAxis];
      if (pos[stepAxis] == out[stepAxis]) {
        break;
      }

      nextCrossingT[stepAxis] += deltaT[stepAxis];
    }

    return hitSomething;
  }

  bool intersectP(Ray ray) {
    // Check ray against overall grid bounds
    List<double> rayT = [0.0];

    if (bounds.inside(ray.pointAt(ray.minDistance))) {
      rayT[0] = ray.minDistance;
    } else if (!bounds.intersectP(ray, rayT)) {
      Stats.GRID_RAY_MISSED_BOUNDS();
      return false;
    }

    Point gridIntersect = ray.pointAt(rayT[0]);

    // Set up 3D DDA for ray
    List<double> nextCrossingT = [0.0, 0.0, 0.0];
    List<double> deltaT = [0.0, 0.0, 0.0];
    List<int> step = [0, 0, 0];
    List<int> out = [0, 0, 0];
    List<int> pos = [0, 0, 0];

    for (int axis = 0; axis < 3; ++axis) {
      // Compute current voxel for axis
      pos[axis] = posToVoxel(gridIntersect, axis);

      if (ray.direction[axis] >= 0.0) {
        // Handle ray with positive direction for voxel stepping
        nextCrossingT[axis] = rayT[0] +
                 (voxelToPos(pos[axis] + 1, axis) - gridIntersect[axis]) /
                 ray.direction[axis];
        deltaT[axis] = width[axis] / ray.direction[axis];
        step[axis] = 1;
        out[axis] = nVoxels[axis];
      } else {
        // Handle ray with negative direction for voxel stepping
        nextCrossingT[axis] = rayT[0] +
                 (voxelToPos(pos[axis], axis) - gridIntersect[axis]) /
                 ray.direction[axis];
        deltaT[axis] = -width[axis] / ray.direction[axis];
        step[axis] = -1;
        out[axis] = -1;
      }
    }

    // Walk grid for shadow ray
    while (true) {
      // Check for intersection in current voxel and advance to next
      int o = offset(pos[0], pos[1], pos[2]);
      _Voxel voxel = voxels[o];

      Stats.GRID_RAY_TRAVERSED_VOXEL(pos, voxel != null ? voxel.size() : 0);
      if (voxel != null && voxel.intersectP(ray)) {
        return true;
      }

      // Advance to next voxel

      // Find _stepAxis_ for stepping to next voxel
      int bits = ((nextCrossingT[0] < nextCrossingT[1]) ? 4 : 0) +
                 ((nextCrossingT[0] < nextCrossingT[2]) ? 2 : 0) +
                 ((nextCrossingT[1] < nextCrossingT[2]) ? 1 : 0);
      const List<int> cmpToAxis = const [ 2, 1, 2, 1, 2, 2, 0, 0 ];
      int stepAxis = cmpToAxis[bits];

      if (ray.maxDistance < nextCrossingT[stepAxis]) {
        break;
      }

      pos[stepAxis] += step[stepAxis];
      if (pos[stepAxis] == out[stepAxis]) {
        break;
      }
      nextCrossingT[stepAxis] += deltaT[stepAxis];
    }

    return false;
  }

  int posToVoxel(Point P, int axis) {
    int v = ((P[axis] - bounds.pMin[axis]) * invWidth[axis]).toInt();
    return v.clamp(0, nVoxels[axis] - 1);
  }

  double voxelToPos(int p, int axis) {
    return bounds.pMin[axis] + p * width[axis];
  }

  int offset(int x, int y, int z) {
    return z * nVoxels[0] * nVoxels[1] + y * nVoxels[0] + x;
  }

  static GridAccel Create(List<Primitive> prims, ParamSet ps) {
    bool refineImmediately = ps.findOneBool('refineimmediately', true);
    return new GridAccel(prims, refineImmediately);
  }

  List<Primitive> primitives = [];
  List<int> nVoxels = [0, 0, 0];
  BBox bounds = new BBox();
  Vector width = new Vector();
  Vector invWidth = new Vector();
  List<_Voxel> voxels;
}

class _Voxel {
  _Voxel(Primitive prim) {
    allCanIntersect = false;
    primitives.add(prim);
  }

  int size() {
    return primitives.length;
  }

  void addPrimitive(Primitive prim) {
    primitives.add(prim);
  }

  bool intersect(Ray ray, Intersection isect) {
    int si = 0;
    int ei = primitives.length;

    // Refine primitives in voxel if needed
    if (!allCanIntersect) {
      for (int i = si; i < ei; ++i) {
        Primitive prim = primitives[i];
        // Refine primitive _prim_ if it's not intersectable
        if (!prim.canIntersect()) {
          List<Primitive> p = [];
          prim.fullyRefine(p);
          assert(p.length > 0);
          if (p.length == 1) {
            primitives[i] = p[0];
          } else {
            primitives[i] = new GridAccel(p, false);
          }
        }
      }
      allCanIntersect = true;
    }

    // Loop over primitives in voxel and find intersections
    bool hitSomething = false;
    for (int i = si; i < ei; ++i) {
      Primitive prim = primitives[i];
      Stats.GRID_RAY_PRIMITIVE_INTERSECTION_TEST(prim);

      if (prim.intersect(ray, isect)) {
        Stats.GRID_RAY_PRIMITIVE_HIT(prim);
        hitSomething = true;
      }
    }

    return hitSomething;
  }

  bool intersectP(Ray ray) {
    Stats.GRID_INTERSECTIONP_TEST(this, ray);
    // Refine primitives in voxel if needed
    if (!allCanIntersect) {
      for (int i = 0; i < primitives.length; ++i) {
        Primitive prim = primitives[i];
        // Refine primitive _prim_ if it's not intersectable
        if (!prim.canIntersect()) {
          List<Primitive> p = [];
          prim.fullyRefine(p);
          assert(p.length > 0);
          if (p.length == 1) {
            primitives[i] = p[0];
          } else {
            primitives[i] = new GridAccel(p, false);
          }
        }
      }

      allCanIntersect = true;
    }

    for (int i = 0; i < primitives.length; ++i) {
      Primitive prim = primitives[i];
      Stats.GRID_RAY_PRIMITIVE_INTERSECTIONP_TEST(prim);
      if (prim.intersectP(ray)) {
        Stats.GRID_RAY_PRIMITIVE_HIT(prim);
        return true;
      }
    }

    return false;
  }

  List<Primitive> primitives = [];
  bool allCanIntersect;
}
