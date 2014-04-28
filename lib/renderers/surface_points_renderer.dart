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
part of renderers;

class SurfacePointsRenderer extends Renderer {
  SurfacePointsRenderer(this.minDist, this.pCamera, this.time, this.filename);

  Future<OutputImage> render(Scene scene) {
    Stopwatch timer = new Stopwatch()..start();
    LogInfo('Starting SurfacePointsRenderer. This may take a few minutes...');

    // Declare shared variables for Poisson point generation
    BBox octBounds = scene.worldBound;
    octBounds.expand(0.001 * Math.pow(octBounds.volume(), 1.0 / 3.0));
    Octree pointOctree = new Octree(octBounds);

    // Create scene bounding sphere to catch rays that leave the scene
    Point sceneCenter = new Point();
    double sceneRadius = scene.worldBound.boundingSphere(sceneCenter);

    Transform objectToWorld = Transform.Translate(sceneCenter - Point.ZERO);
    Transform worldToObject = Transform.Inverse(objectToWorld);

    Shape sph = new Sphere(objectToWorld, worldToObject,
                           true, sceneRadius, -sceneRadius,
                           sceneRadius, 360.0);

    GeometricPrimitive sphere = new GeometricPrimitive(sph, null, null);

    _SurfacePointStats stats = new _SurfacePointStats();
    int maxFails = 2000;

    if (RenderOverrides.QuickRender()) {
      maxFails = Math.max(10, maxFails ~/ 10);
    }

    // Launch tasks to trace rays to find Poisson points
    Stats.SUBSURFACE_STARTED_RAYS_FOR_POINTS();

    int nTasks = 1;
    for (int i = 0; i < nTasks; ++i) {
      _SurfacePointTask task = new _SurfacePointTask(scene, pCamera, time, i,
                                                     minDist, maxFails, stats,
                                                     sphere, pointOctree,
                                                     points);
      task.run();
    }

    LogInfo('FINISHED SurfacePointsRenderer, ${points.length} points : '
            '${timer.elapsed}');

    Stats.SUBSURFACE_FINISHED_RAYS_FOR_POINTS(stats.totalRaysTraced,
                                              stats.numPointsAdded);

    if (filename.isNotEmpty) {
      Float32List data = new Float32List(points.length * 8);
      int di = 0;
      for (int i = 0; i < points.length; ++i) {
        SurfacePoint sp = points[i];
        data[di++] = sp.p.x;
        data[di++] = sp.p.y;
        data[di++] = sp.p.z;
        data[di++] = sp.n.x;
        data[di++] = sp.n.y;
        data[di++] = sp.n.z;
        data[di++] = sp.area;
        data[di++] = sp.rayEpsilon;
      }

      ResourceManager.WriteFile(filename, data);
    }

    // This renderer does not generate an image.
    Completer<OutputImage> c = new Completer<OutputImage>();
    c.complete(null);
    return c.future;
  }

  Spectrum Li(Scene scene, RayDifferential ray, Sample sample, RNG rng,
              [Intersection isect, Spectrum T]) {
    return new Spectrum(0.0);
  }

  Spectrum transmittance(Scene scene, RayDifferential ray,
                         Sample sample, RNG rng) {
    return new Spectrum(0.0);
  }

  double minDist;
  double time;
  Point pCamera;
  String filename;
  List<SurfacePoint> points = [];

  static void FindPoissonPointDistribution(Point pCamera, double time,
                                             double minDist, Scene scene,
                                             List<SurfacePoint> points) {
    SurfacePointsRenderer sp = new SurfacePointsRenderer(minDist, pCamera,
                                                         time, '');
    sp.render(scene);
    points.clear();
    points.insertAll(0, sp.points);
  }

  static SurfacePointsRenderer Create(ParamSet params, Point pCamera,
                                      double time) {
    double minDist = params.findOneFloat('minsampledistance', 0.25);
    String filename = params.findOneFilename('filename', '');
    if (RenderOverrides.QuickRender()) {
      minDist *= 4.0;
    }
    return new SurfacePointsRenderer(minDist, pCamera, time, filename);
  }
}

class _SurfacePointStats {
  int repeatedFails = 0;
  int maxRepeatedFails = 0;
  int totalPathsTraced = 0;
  int totalRaysTraced = 0;
  int numPointsAdded = 0;
}

class _PoissonCheck {
  _PoissonCheck(double md, this.p) :
    maxDist2 = md * md,
    failed = false;

  bool call(SurfacePoint sp) {
    if (Vector.DistanceSquared(sp.p, p) < maxDist2) {
      failed = true;
      return false;
    }
    return true;
  }

  double maxDist2;
  bool failed;
  Point p;
}


class _SurfacePointTask {
  _SurfacePointTask(this.scene, this.origin, this.time, this.taskNum,
                    this.minSampleDist, this.maxFails, this.stats,
                    this.sphere, this.octree, this.surfacePoints);

  void run() {
    RNG rng = new RNG(37 * taskNum);
    List<SurfacePoint> candidates = [];

    Intersection isect = new Intersection();

    while (true) {
      int pathsTraced;
      int raysTraced = 0;

      for (pathsTraced = 0; pathsTraced < 20000; ++pathsTraced) {
        // Follow ray path and attempt to deposit candidate sample points
        Vector dir = UniformSampleSphere(rng.randomFloat(), rng.randomFloat());
        Ray ray = new Ray(origin, dir, 0.0, INFINITY, time);

        while (ray.depth < 30) {
          // Find ray intersection with scene geometry or bounding sphere
          ++raysTraced;
          bool hitOnSphere = false;

          if (!scene.intersect(ray, isect)) {
            if (!sphere.intersect(ray, isect)) {
              break;
            }
            hitOnSphere = true;
          }

          DifferentialGeometry hitGeometry = isect.dg;
          hitGeometry.nn = Normal.FaceForward(hitGeometry.nn, -ray.direction);

          // Store candidate sample point at ray intersection if appropriate
          if (!hitOnSphere && ray.depth >= 3 &&
              isect.getBSSRDF(new RayDifferential.fromRay(ray)) != null) {
            double area = Math.PI * (minSampleDist / 2.0) *
                          (minSampleDist / 2.0);
            candidates.add(new SurfacePoint(hitGeometry.p, hitGeometry.nn,
                                            area, isect.rayEpsilon));
          }

          // Generate random ray from intersection point
          Vector dir = UniformSampleSphere(rng.randomFloat(), rng.randomFloat());
          dir = Vector.FaceForward(dir, hitGeometry.nn);
          ray = new Ray.withParent(hitGeometry.p, dir, ray, isect.rayEpsilon);
        }
      }

      // Make first pass through candidate points with reader lock
      List<bool> candidateRejected = [];
      for (int i = 0; i < candidates.length; ++i) {
        _PoissonCheck check = new _PoissonCheck(minSampleDist, candidates[i].p);
        octree.lookup(candidates[i].p, check);
        candidateRejected.add(check.failed);
      }

      // Make second pass through points with writer lock and update octree
      if (stats.repeatedFails >= maxFails) {
        return;
      }

      stats.totalPathsTraced += pathsTraced;
      stats.totalRaysTraced += raysTraced;
      int oldMaxRepeatedFails = stats.maxRepeatedFails;

      for (int i = 0; i < candidates.length; ++i) {
        if (candidateRejected[i]) {
          // Update for rejected candidate point
          ++stats.repeatedFails;
          stats.maxRepeatedFails = Math.max(stats.maxRepeatedFails,
                                            stats.repeatedFails);
          if (stats.repeatedFails >= maxFails) {
            return;
          }
        } else {
          // Recheck candidate point and possibly add to octree
          SurfacePoint sp = candidates[i];
          _PoissonCheck check = new _PoissonCheck(minSampleDist, sp.p);
          octree.lookup(sp.p, check);
          if (check.failed) {
            // Update for rejected candidate point
            ++stats.repeatedFails;
            stats.maxRepeatedFails = Math.max(stats.maxRepeatedFails,
                                              stats.repeatedFails);
            if (stats.repeatedFails >= maxFails) {
              return;
            }
          } else {
            ++stats.numPointsAdded;
            stats.repeatedFails = 0;
            Vector delta = new Vector(minSampleDist, minSampleDist,
                                      minSampleDist);
            octree.add(sp, new BBox(sp.p - delta, sp.p + delta));
            Stats.SUBSURFACE_ADDED_POINT_TO_OCTREE(sp, minSampleDist);
            surfacePoints.add(sp);
          }
        }
      }

      /*LogDebug('SurfacePointsRenderer: Tracing Paths: $raysTraced : '
               '${stats.repeatedFails} / ${maxFails} : '
               '${surfacePoints.length}');*/

      // Stop following paths if not finding new points
      if (stats.totalPathsTraced > 50000 && stats.numPointsAdded == 0) {
        LogWarning('There don\'t seem to be any objects with BSSRDFs '
                   'in this scene. Giving up.');
        return;
      }

      candidates.clear();
    }
  }

  int taskNum;
  Scene scene;
  Point origin;
  double time;
  double minSampleDist;
  int maxFails;

  _SurfacePointStats stats;
  GeometricPrimitive sphere;
  Octree octree;
  List<SurfacePoint> surfacePoints = [];
}
