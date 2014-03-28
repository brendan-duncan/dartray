part of core;

/**
 * Collect statistics about the rendering process.
 */
class Stats {
  static String getString() {
    return '';
  }

  static void CREATED_SHAPE(Shape shape) {
  }

  static void CREATED_TRIANGLE(tri) {
  }

  static void STARTED_GENERATING_CAMERA_RAY(CameraSample s) {
  }

  static void KDTREE_CREATED_INTERIOR_NODE(int axis, double split) {
  }

  static void KDTREE_CREATED_LEAF(int nprims, int depth) {
  }

  static void RAY_TRIANGLE_INTERSECTION_TEST(Ray ray, tri) {
  }

  static void RAY_TRIANGLE_INTERSECTIONP_TEST(Ray ray, tri) {
  }

  static void RAY_TRIANGLE_INTERSECTION_HIT(Ray ray, double t) {
  }

  static void RAY_TRIANGLE_INTERSECTIONP_HIT(Ray ray, double t) {
  }

  static void FINISHED_RAY_INTERSECTION(Ray ray, Intersection isect, bool hit) {
  }

  static void FINISHED_RAY_INTERSECTIONP(Ray ray, bool hit) {
  }

  static void STARTED_SPECULAR_REFLECTION_RAY(RayDifferential ray) {
  }

  static void STARTED_SPECULAR_REFRACTION_RAY(RayDifferential ray) {
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

  static void LOADED_IMAGE_MAP(arg0, arg1, arg2, arg3, arg4) {
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
}
