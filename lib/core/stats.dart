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
 * Collect statistics about the rendering process.
 */
class Stats {
  static String getString() {
    String s = '';
    for (StatTracker t in trackers) {
      s += '$t\n';
    }
    return s;
  }

  static void CREATED_SHAPE(Shape shape) {
    shapesMade++;
  }

  static void CREATED_TRIANGLE(tri) {
    trianglesMade++;
  }

  static void STARTED_GENERATING_CAMERA_RAY(CameraSample s) {
    cameraRays++;
  }

  static void KDTREE_CREATED_INTERIOR_NODE(int axis, double split) {
    kdTreeInteriorNodes++;
  }

  static void KDTREE_CREATED_LEAF(int nprims, int depth) {
    kdTreeLeafNodes++;
    kdTreeMaxPrims.max(nprims);
    kdTreeMaxDepth.max(depth);
  }

  static void RAY_TRIANGLE_INTERSECTION_TEST(Ray ray, tri) {
    rayTriIntersections.add(0, 1);
  }

  static void RAY_TRIANGLE_INTERSECTIONP_TEST(Ray ray, tri) {
    rayTriIntersectionPs.add(0, 1);
  }

  static void RAY_TRIANGLE_INTERSECTION_HIT(Ray ray, double t) {
    rayTriIntersections.add(1, 0);
  }

  static void RAY_TRIANGLE_INTERSECTIONP_HIT(Ray ray, double t) {
    rayTriIntersectionPs.add(1, 0);
  }

  static void FINISHED_RAY_INTERSECTION(Ray ray, Intersection isect, bool hit) {
    nonShadowRays++;
  }

  static void FINISHED_RAY_INTERSECTIONP(Ray ray, bool hit) {
    shadowRays++;
  }

  static void STARTED_SPECULAR_REFLECTION_RAY(RayDifferential ray) {
    specularReflectionRays++;
  }

  static void STARTED_SPECULAR_REFRACTION_RAY(RayDifferential ray) {
    specularRefractionRays++;
  }

  static void ACCESSED_TEXEL(arg0, arg1, arg2, arg3) {
  }

  static void ALLOCATED_CACHED_TRANSFORM() {
  }

  static void FOUND_CACHED_TRANSFORM() {
  }

  static void ATOMIC_MEMORY_OP() {
  }

  static void BVH_STARTED_CONSTRUCTION(arg0, arg1) {
  }

  static void BVH_FINISHED_CONSTRUCTION(arg0) {
  }

  static void BVH_INTERSECTION_STARTED(arg0, arg1) {
  }

  static void BVH_INTERSECTION_TRAVERSED_INTERIOR_NODE(arg0) {
  }

  static void BVH_INTERSECTION_TRAVERSED_LEAF_NODE(arg0) {
  }

  static void BVH_INTERSECTION_PRIMITIVE_TEST(arg0) {
  }

  static void BVH_INTERSECTION_PRIMITIVE_HIT(arg0) {
  }

  static void BVH_INTERSECTION_PRIMITIVE_MISSED(arg0) {
  }

  static void BVH_INTERSECTION_FINISHED() {
  }

  static void BVH_INTERSECTIONP_STARTED(arg0, arg1) {
  }

  static void BVH_INTERSECTIONP_TRAVERSED_INTERIOR_NODE(arg0) {
  }

  static void BVH_INTERSECTIONP_TRAVERSED_LEAF_NODE(arg0) {
  }

  static void BVH_INTERSECTIONP_PRIMITIVE_TEST(arg0) {
  }

  static void BVH_INTERSECTIONP_PRIMITIVE_HIT(arg0) {
  }

  static void BVH_INTERSECTIONP_PRIMITIVE_MISSED(arg0) {
  }

  static void BVH_INTERSECTIONP_FINISHED() {
  }

  static void FINISHED_PARSING() {
  }

  static void FINISHED_PREPROCESSING() {
  }

  static void FINISHED_RENDERING() {
  }

  static void FINISHED_RENDERTASK(arg0) {
  }

  static void FINISHED_TASK(arg0) {
  }

  static void FINISHED_ADDING_IMAGE_SAMPLE() {
  }

  static void FINISHED_CAMERA_RAY_INTEGRATION(arg0, arg1, arg2) {
  }

  static void FINISHED_EWA_TEXTURE_LOOKUP() {
  }

  static void FINISHED_GENERATING_CAMERA_RAY(arg0, arg1, arg2) {
  }

  static void FINISHED_BSDF_SHADING(arg0, arg1) {
  }

  static void FINISHED_BSSRDF_SHADING(arg0, arg1) {
  }

  static void FINISHED_SPECULAR_REFLECTION_RAY(arg0) {
  }

  static void FINISHED_SPECULAR_REFRACTION_RAY(arg0) {
  }

  static void FINISHED_TRILINEAR_TEXTURE_LOOKUP() {
  }

  static void GRID_BOUNDS_AND_RESOLUTION(arg0, arg1) {
  }

  static void GRID_FINISHED_CONSTRUCTION(arg0) {
  }

  static void GRID_INTERSECTIONP_TEST(arg0, arg1) {
  }

  static void GRID_INTERSECTION_TEST(arg0, arg1) {
  }

  static void GRID_RAY_MISSED_BOUNDS() {
  }

  static void GRID_RAY_PRIMITIVE_HIT(arg0) {
  }

  static void GRID_RAY_PRIMITIVE_INTERSECTIONP_TEST(arg0) {
  }

  static void GRID_RAY_PRIMITIVE_INTERSECTION_TEST(arg0) {
  }

  static void GRID_RAY_TRAVERSED_VOXEL(arg0, arg1) {
  }

  static void GRID_STARTED_CONSTRUCTION(arg0, arg1) {
  }

  static void GRID_VOXELIZED_PRIMITIVE(arg0, arg1) {
  }

  static void IRRADIANCE_CACHE_ADDED_NEW_SAMPLE(arg0, arg1, arg2, arg3, arg4, arg5) {
  }

  static void IRRADIANCE_CACHE_CHECKED_SAMPLE(arg0, arg1, arg2) {
  }

  static void IRRADIANCE_CACHE_FINISHED_COMPUTING_IRRADIANCE(arg0, arg1) {
  }

  static void IRRADIANCE_CACHE_FINISHED_INTERPOLATION(arg0, arg1, arg2, arg3) {
  }

  static void IRRADIANCE_CACHE_FINISHED_RAY(arg0, arg1, arg2) {
  }

  static void IRRADIANCE_CACHE_STARTED_COMPUTING_IRRADIANCE(arg0, arg1) {
  }

  static void IRRADIANCE_CACHE_STARTED_INTERPOLATION(arg0, arg1) {
  }

  static void IRRADIANCE_CACHE_STARTED_RAY(arg0) {
  }

  static void KDTREE_FINISHED_CONSTRUCTION(arg0) {
  }

  static void KDTREE_INTERSECTIONP_PRIMITIVE_TEST(arg0) {
  }

  static void KDTREE_INTERSECTION_PRIMITIVE_TEST(arg0) {
  }

  static void KDTREE_INTERSECTIONP_HIT(arg0) {
  }

  static void KDTREE_INTERSECTIONP_MISSED() {
  }

  static void KDTREE_INTERSECTIONP_TEST(arg0, arg1) {
  }

  static void KDTREE_INTERSECTION_FINISHED() {
  }

  static void KDTREE_INTERSECTION_HIT(arg0) {
  }

  static void KDTREE_INTERSECTION_TEST(arg0, arg1) {
  }

  static void KDTREE_RAY_MISSED_BOUNDS() {
  }

  static void KDTREE_STARTED_CONSTRUCTION(arg0, arg1) {
  }

  static void KDTREE_INTERSECTION_TRAVERSED_INTERIOR_NODE(arg0) {
  }

  static void KDTREE_INTERSECTION_TRAVERSED_LEAF_NODE(arg0, arg1) {
  }

  static void KDTREE_INTERSECTIONP_TRAVERSED_INTERIOR_NODE(arg0) {
  }

  static void KDTREE_INTERSECTIONP_TRAVERSED_LEAF_NODE(arg0, arg1) {
  }

  static void LOADED_IMAGE_MAP(arg0, arg1, arg2, arg3) {
  }

  static void MIPMAP_EWA_FILTER(arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) {
  }

  static void MIPMAP_TRILINEAR_FILTER(arg0, arg1, arg2, arg3, arg4, arg5) {
  }

  static void MLT_ACCEPTED_MUTATION(arg0, arg1, arg2) {
  }

  static void MLT_REJECTED_MUTATION(arg0, arg1, arg2) {
  }

  static void MLT_STARTED_MLT_TASK(arg0) {
  }

  static void MLT_FINISHED_MLT_TASK(arg0) {
  }

  static void MLT_STARTED_RENDERING() {
  }

  static void MLT_FINISHED_RENDERING() {
  }

  static void MLT_STARTED_DIRECTLIGHTING() {
  }

  static void MLT_FINISHED_DIRECTLIGHTING() {
  }

  static void MLT_STARTED_BOOTSTRAPPING(count) {
  }

  static void MLT_FINISHED_BOOTSTRAPPING(b) {
  }

  static void MLT_STARTED_MUTATION() {
  }

  static void MLT_FINISHED_MUTATION() {
  }

  static void MLT_STARTED_SAMPLE_SPLAT() {
  }

  static void MLT_FINISHED_SAMPLE_SPLAT() {
  }

  static void MLT_STARTED_GENERATE_PATH() {
  }

  static void MLT_FINISHED_GENERATE_PATH() {
  }

  static void MLT_STARTED_LPATH() {
  }

  static void MLT_FINISHED_LPATH() {
  }

  static void MLT_STARTED_LBIDIR() {
  }

  static void MLT_FINISHED_LBIDIR() {
  }

  static void MLT_STARTED_TASK_INIT() {
  }

  static void MLT_FINISHED_TASK_INIT() {
  }

  static void MLT_STARTED_SAMPLE_LIGHT_FOR_BIDIR() {
  }

  static void MLT_FINISHED_SAMPLE_LIGHT_FOR_BIDIR() {
  }

  static void MLT_STARTED_DISPLAY_UPDATE() {
  }

  static void MLT_FINISHED_DISPLAY_UPDATE() {
  }

  static void MLT_STARTED_ESTIMATE_DIRECT() {
  }

  static void MLT_FINISHED_ESTIMATE_DIRECT() {
  }

  static void PHOTON_MAP_DEPOSITED_CAUSTIC_PHOTON(arg0, arg1, arg2) {
  }

  static void PHOTON_MAP_DEPOSITED_DIRECT_PHOTON(arg0, arg1, arg2) {
  }

  static void PHOTON_MAP_DEPOSITED_INDIRECT_PHOTON(arg0, arg1, arg2) {
  }

  static void PHOTON_MAP_FINISHED_GATHER_RAY(arg0) {
  }

  static void PHOTON_MAP_FINISHED_LOOKUP(arg0, arg1, arg2, arg3) {
  }

  static void PHOTON_MAP_FINISHED_RAY_PATH(arg0, arg1) {
  }

  static void PHOTON_MAP_STARTED_GATHER_RAY(arg0) {
  }

  static void PHOTON_MAP_STARTED_LOOKUP(arg0) {
  }

  static void PHOTON_MAP_STARTED_RAY_PATH(arg0, arg1) {
  }

  static void SAMPLE_OUTSIDE_IMAGE_EXTENT(arg0) {
  }

  static void STARTED_ADDING_IMAGE_SAMPLE(arg0, arg1, arg2, arg3) {
  }

  static void STARTED_CAMERA_RAY_INTEGRATION(arg0, arg1) {
  }

  static void STARTED_EWA_TEXTURE_LOOKUP(arg0, arg1) {
  }

  static void STARTED_PARSING() {
  }

  static void STARTED_PREPROCESSING() {
  }

  static void STARTED_RAY_INTERSECTION(arg0) {
  }

  static void STARTED_RAY_INTERSECTIONP(arg0) {
  }

  static void STARTED_RENDERING() {
  }

  static void STARTED_RENDERTASK(arg0) {
  }

  static void STARTED_BSDF_SHADING(arg0) {
  }

  static void STARTED_BSSRDF_SHADING(arg0) {
  }

  static void STARTED_TASK(arg0) {
  }

  static void STARTED_TRILINEAR_TEXTURE_LOOKUP(arg0, arg1) {
  }

  static void SUBSURFACE_ADDED_INTERIOR_CONTRIBUTION(arg0) {
  }

  static void SUBSURFACE_ADDED_POINT_CONTRIBUTION(arg0) {
  }

  static void SUBSURFACE_ADDED_POINT_TO_OCTREE(arg0, arg1) {
  }

  static void SUBSURFACE_COMPUTED_IRRADIANCE_AT_POINT(arg0, arg1) {
  }

  static void SUBSURFACE_FINISHED_COMPUTING_IRRADIANCE_VALUES() {
  }

  static void SUBSURFACE_FINISHED_OCTREE_LOOKUP() {
  }

  static void SUBSURFACE_FINISHED_RAYS_FOR_POINTS(arg0, arg1) {
  }

  static void SUBSURFACE_STARTED_COMPUTING_IRRADIANCE_VALUES() {
  }

  static void SUBSURFACE_STARTED_OCTREE_LOOKUP(arg0) {
  }

  static void SUBSURFACE_STARTED_RAYS_FOR_POINTS() {
  }

  static void SUPERSAMPLE_PIXEL_NO(arg0, arg1) {
  }

  static void SUPERSAMPLE_PIXEL_YES(arg0, arg1) {
  }

  static void RNG_STARTED_RANDOM_FLOAT() {
  }

  static void RNG_FINISHED_RANDOM_FLOAT() {
  }

  static void RNG_FINISHED_TABLEGEN() {
  }

  static void RNG_STARTED_TABLEGEN() {
  }

  static void STARTED_BSDF_EVAL() {
  }

  static void FINISHED_BSDF_EVAL() {
  }

  static void STARTED_BSDF_SAMPLE() {
  }

  static void FINISHED_BSDF_SAMPLE() {
  }

  static void STARTED_BSDF_PDF() {
  }

  static void FINISHED_BSDF_PDF() {
  }

  static void AREA_LIGHT_STARTED_SAMPLE() {
  }

  static void AREA_LIGHT_FINISHED_SAMPLE() {
  }

  static void INFINITE_LIGHT_STARTED_SAMPLE() {
  }

  static void INFINITE_LIGHT_FINISHED_SAMPLE() {
  }

  static void INFINITE_LIGHT_STARTED_PDF() {
  }

  static void INFINITE_LIGHT_FINISHED_PDF() {
  }

  static List<StatTracker> trackers = [];

  static StatsCounter shapesMade = new StatsCounter('Shapes', 'Total Shapes Created');
  static StatsCounter trianglesMade = new StatsCounter('Shapes', 'Total Triangles Created');
  static StatsCounter cameraRays = new StatsCounter('Rays', 'Camera Rays Traced');
  static StatsCounter specularReflectionRays = new StatsCounter('Rays', 'Specular Reflection Rays Traced');
  static StatsCounter specularRefractionRays = new StatsCounter('Rays', 'Specular Refraction Rays Traced');
  static StatsCounter shadowRays = new StatsCounter('Rays', 'Shadow Rays Traced');
  static StatsCounter nonShadowRays = new StatsCounter('Rays', 'Total Non-Shadow Rays Traced');
  static StatsCounter kdTreeInteriorNodes = new StatsCounter('Kd-Tree', 'Interior Nodes Created');
  static StatsCounter kdTreeLeafNodes = new StatsCounter('Kd-Tree', 'Interior Nodes Created');
  static StatsCounter kdTreeMaxPrims = new StatsCounter('Kd-Tree', 'Maximum Primitives in Leaf');
  static StatsCounter kdTreeMaxDepth = new StatsCounter('Kd-Tree', 'Maximum Depth of Leaf Nodes');
  static StatsPercentage rayTriIntersections = new StatsPercentage('Intersections', 'Ray/Triangle Intersection Hits');
  static StatsPercentage rayTriIntersectionPs = new StatsPercentage('Intersections', 'Ray/Triangle IntersectionP Hits');
}

class StatTracker {
  String category;
  String name;

  StatTracker(this.category, this.name);
}

class StatsCounter extends StatTracker {
  StatsCounter(String category, String name) :
    super(category, name),
    count = 0 {
    Stats.trackers.add(this);
  }

  StatsCounter operator +(int n) {
    count += n;
    return this;
  }

  void max(int newval) {
    count = Math.max(count, newval);
  }

  void min(int newval) {
    count = Math.min(count, newval);
  }

  String toString() => '$category | $name: $count';

  int count;
}

class StatsPercentage extends StatTracker {
  StatsPercentage(String category, String name) :
    super(category, name),
    na = 0,
    nb = 0 {
    Stats.trackers.add(this);
  }

  void add(int a, int b) {
    na += a;
    nb += b;
  }

  String toString() => '$category | $name: $na / $nb (${(na / nb) * 100}';

  int na;
  int nb;
}
