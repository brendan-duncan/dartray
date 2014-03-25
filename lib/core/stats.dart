part of core;

/**
 * Collect statistics about the rendering process.
 */
class Stats {
  static int rayIntersection = 0;
  static int rayIntersectionP = 0;

  static void StartedRayIntersection(Ray r) {
  }

  static void FinishedRayIntersection(Ray r) {
    rayIntersection++;
  }

  static void StartedRayIntersectionP(Ray r) {
  }

  static void FinishedRayIntersectionP(Ray r) {
    rayIntersectionP++;
  }
}
