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

class CreateProbesRenderer extends Renderer {
  CreateProbesRenderer(this.surfaceIntegrator, this.volumeIntegrator,
                       this.camera, this.lmax, this.probeSpacing, this.bbox,
                       this.nIndirSamples, this.includeDirectInProbes,
                       this.includeIndirectInProbes, this.time, this.filename);

  Future<OutputImage> render(Scene scene) {
    LogInfo('Starting CreateProbesRenderer');
    // Compute scene bounds and initialize probe integrators
    if (bbox.pMin.x > bbox.pMax.x) {
      bbox = scene.worldBound;
    }

    surfaceIntegrator.preprocess(scene, camera, this);
    volumeIntegrator.preprocess(scene, camera, this);

    Sample origSample = new Sample(null, surfaceIntegrator, volumeIntegrator,
                                   scene);

    // Compute sampling rate in each dimension
    Vector delta = bbox.pMax - bbox.pMin;
    List<int> nProbes = [0, 0, 0];
    for (int i = 0; i < 3; ++i) {
      nProbes[i] = Math.max(1, (delta[i] / probeSpacing).ceil());
    }

    // Allocate SH coefficient vector pointers for sample points
    int count = nProbes[0] * nProbes[1] * nProbes[2];
    List<List<Spectrum>> c_in = new List<List<Spectrum>>(count);
    for (int i = 0; i < count; ++i) {
      c_in[i] = new List<Spectrum>(SphericalHarmonics.Terms(lmax));
    }

    // Compute random points on surfaces of scene

    // Create scene bounding sphere to catch rays that leave the scene
    Point sceneCenter = new Point();
    double sceneRadius = scene.worldBound.boundingSphere(sceneCenter);
    Transform objectToWorld = Transform.Translate(sceneCenter - Point.ZERO);
    Transform worldToObject = Transform.Inverse(objectToWorld);
    Shape sph = new Sphere(objectToWorld, worldToObject, true, sceneRadius,
                           -sceneRadius, sceneRadius, 360.0);
    GeometricPrimitive sphere = new GeometricPrimitive(sph, null, null);
    const int nPoints = 32768;
    const int maxDepth = 32;
    List<Point> surfacePoints = [];

    Point pCamera = camera.cameraToWorld.transformPoint(camera.shutterOpen,
                                                        Point.ZERO);
    surfacePoints.add(pCamera);

    RNG rng = new RNG();

    while (surfacePoints.length < nPoints) {
      // Generate random path from camera and deposit surface points
      Point pray = pCamera;
      Vector dir = UniformSampleSphere(rng.randomFloat(), rng.randomFloat());
      double rayEpsilon = 0.0;
      for (int i = 0; i < maxDepth; ++i) {
        Ray ray = new Ray(pray, dir, rayEpsilon, INFINITY, time);

        Intersection isect = new Intersection();
        if (!scene.intersect(ray, isect) && !sphere.intersect(ray, isect)) {
          break;
        }

        surfacePoints.add(ray.pointAt(ray.maxDistance));

        DifferentialGeometry hitGeometry = isect.dg;
        pray = isect.dg.p;
        rayEpsilon = isect.rayEpsilon;
        hitGeometry.nn = Normal.FaceForward(hitGeometry.nn, -ray.direction);

        dir = UniformSampleSphere(rng.randomFloat(), rng.randomFloat());
        dir = Vector.FaceForward(dir, hitGeometry.nn);
      }
    }

    // Launch tasks to compute radiance probes at sample points
    for (int i = 0; i < count; ++i) {
      _CreateRadProbeTask task = new _CreateRadProbeTask(i, nProbes, time,
                                    bbox, lmax, includeDirectInProbes,
                                    includeIndirectInProbes, nIndirSamples,
                                    origSample, surfacePoints,
                                    scene, this, c_in[i]);
      task.run();
    }

    // Write radiance probe coefficients to file
    int size = 15 + (nProbes[0] * nProbes[1] * nProbes[2]) *
               SphericalHarmonics.Terms(lmax) * Spectrum.NumSamples();

    // Store the generated data as a resource so that subsequent renderers can
    // retrieve it.
    Float32List fp = new Float32List(size);
    int fi = 0;
    fp[fi++] = lmax.toDouble();
    fp[fi++] = includeDirectInProbes ? 1.0 : 0.0;
    fp[fi++] = includeIndirectInProbes ? 1.0 : 0.0;
    fp[fi++] = nProbes[0].toDouble();
    fp[fi++] = nProbes[1].toDouble();
    fp[fi++] = nProbes[2].toDouble();
    fp[fi++] = bbox.pMin.x;
    fp[fi++] = bbox.pMin.y;
    fp[fi++] = bbox.pMin.z;
    fp[fi++] = bbox.pMax.x;
    fp[fi++] = bbox.pMax.y;
    fp[fi++] = bbox.pMax.z;
    for (int i = 0, il = nProbes[0] * nProbes[1] * nProbes[2]; i < il; ++i) {
      for (int j = 0, jl = SphericalHarmonics.Terms(lmax); j < jl; ++j) {
        List<double> s = c_in[i][j].toList();
        for (int k = 0; k < s.length; ++k) {
          fp[fi++] = s[k];
        }
      }
    }

    ResourceManager.WriteFile(filename, fp);

    // This renderer does not generate an image.
    Completer<OutputImage> c = new Completer<OutputImage>();
    c.complete(null);
    return c.future;
  }

  Spectrum Li(Scene scene, RayDifferential ray, Sample sample, RNG rng,
              [Intersection isect, Spectrum T]) {
    assert(ray.time == sample.time);

    if (T == null) {
      T = new Spectrum();
    }

    if (isect == null) {
      isect = new Intersection();
    }

    assert(!ray.hasNaNs());

    Spectrum Lo = new Spectrum(0.0);

    if (scene.intersect(ray, isect)) {
      Lo = surfaceIntegrator.Li(scene, this, ray, isect, sample, rng);
    } else {
      for (int i = 0; i < scene.lights.length; ++i) {
        Lo += scene.lights[i].Le(ray);
      }
    }

    Spectrum Lv = volumeIntegrator.Li(scene, this, ray, sample, rng, T);

    return T * Lo + Lv;
  }

  Spectrum transmittance(Scene scene, RayDifferential ray,
                         Sample sample, RNG rng) {
    return volumeIntegrator.transmittance(scene, this, ray, sample, rng);
  }

  static CreateProbesRenderer Create(Camera camera, SurfaceIntegrator surf,
                                     VolumeIntegrator vol, ParamSet params) {
    bool includeDirect = params.findOneBool('directlighting', true);
    bool includeIndirect = params.findOneBool('indirectlighting', true);
    int lmax = params.findOneInt('lmax', 4);
    int nindir = params.findOneInt('indirectsamples', 512);
    int nbbox;
    BBox bounds;
    List<double> b = params.findFloat('bounds');

    if (b != null) {
      if (b.length != 6) {
        LogWarning('Expecting six values [x0 y0 z0 x1 y1 z1] for bounds');
      } else {
        bounds = new BBox(new Point(b[0], b[1], b[2]),
                          new Point(b[3], b[4], b[5]));
      }
    }

    double probeSpacing = params.findOneFloat('samplespacing', 1.0);
    double time = params.findOneFloat('time', 0.0);
    String filename = params.findOneFilename('filename', 'probes.out');

    return new CreateProbesRenderer(surf, vol, camera, lmax, probeSpacing,
                                    bounds, nindir, includeDirect,
                                    includeIndirect, time, filename);
  }

  SurfaceIntegrator surfaceIntegrator;
  VolumeIntegrator volumeIntegrator;
  Camera camera;
  int lmax;
  int nIndirSamples;
  BBox bbox;
  bool includeDirectInProbes;
  bool includeIndirectInProbes;
  double time;
  double probeSpacing;
  String filename;
}


class _CreateRadProbeTask {
  _CreateRadProbeTask(this.pointNum, List<int> nProbes, this.time, this.bbox,
                      this.lmax, this.includeDirectInProbes,
                      this.includeIndirectInProbes, this.nIndirSamples,
                      this.origSample, this.surfacePoints, this.scene,
                      this.renderer, this.c_in) :
    this.nProbes = new List<int>.from(nProbes);

  void run() {
    // Compute region in which to compute incident radiance probes
    int sx = pointNum % nProbes[0];
    int sy = (pointNum ~/ nProbes[0]) % nProbes[1];
    int sz = pointNum ~/ (nProbes[0] * nProbes[1]);
    assert(sx >= 0 && sx < nProbes[0]);
    assert(sy >= 0 && sy < nProbes[1]);
    assert(sz >= 0 && sz < nProbes[2]);
    double tx0 = sx / nProbes[0];
    double tx1 = (sx + 1) / nProbes[0];
    double ty0 = sy / nProbes[1];
    double ty1 = (sy + 1) / nProbes[1];
    double tz0 = sz / nProbes[2];
    double tz1 = (sz + 1) / nProbes[2];
    BBox b = new BBox(bbox.lerp(tx0, ty0, tz0), bbox.lerp(tx1, ty1, tz1));

    final lmax_terms = SphericalHarmonics.Terms(lmax);

    // Initialize common variables for _CreateRadProbeTask::Run()_
    RNG rng = new RNG(pointNum);
    List<Spectrum> c_probe = new List<Spectrum>(lmax_terms);
    int nFound = 0;
    int lastVisibleOffset = 0;

    for (int i = 0; i < 256; ++i) {
      if (nFound == 32) {
        break;
      }

      // Try to compute radiance probe contribution at _i_th sample point

      // Compute _i_th candidate point _p_ in cell's bounding box
      double dx = RadicalInverse(i + 1, 2);
      double dy = RadicalInverse(i + 1, 3);
      double dz = RadicalInverse(i + 1, 5);
      Point p = b.lerp(dx, dy, dz);

      // Skip point _p_ if not indirectly visible from camera
      if (scene.intersectP(new Ray(surfacePoints[lastVisibleOffset],
                                   p - surfacePoints[lastVisibleOffset],
                                   1e-4, 1.0, time))) {
        int j;
        // See if point is visible to any element of _surfacePoints_
        for (j = 0; j < surfacePoints.length; ++j) {
          if (!scene.intersectP(new Ray(surfacePoints[j], p - surfacePoints[j],
                                        1e-4, 1.0, time))) {
            lastVisibleOffset = j;
            break;
          }
        }
        if (j == surfacePoints.length) {
          continue;
        }
      }

      ++nFound;

      // Compute SH coefficients of incident radiance at point _p_
      if (includeDirectInProbes) {
        for (int i = 0; i < lmax_terms; ++i) {
          c_probe[i] = new Spectrum(0.0);
        }

        SphericalHarmonics.ProjectIncidentDirectRadiance(p, 0.0, time, scene,
                                                         true, lmax, rng,
                                                         c_probe);
        for (int i = 0; i < lmax_terms; ++i) {
          c_in[i] += c_probe[i];
        }
      }

      if (includeIndirectInProbes) {
        for (int i = 0; i < lmax_terms; ++i) {
          c_probe[i] = new Spectrum(0.0);
        }

        SphericalHarmonics.ProjectIncidentIndirectRadiance(p, 0.0, time,
                                                           renderer,
                                                           origSample, scene,
                                                           lmax, rng,
                                                           nIndirSamples,
                                                           c_probe);
        for (int i = 0; i < lmax_terms; ++i) {
          c_in[i] += c_probe[i];
        }
      }
    }

    // Compute final average value for probe and cleanup
    if (nFound > 0) {
      for (int i = 0; i < lmax_terms; ++i) {
        c_in[i] /= nFound;
      }
    }
  }

  int pointNum;
  List<int> nProbes;
  BBox bbox;
  int lmax;
  int nIndirSamples;
  double time;
  bool includeDirectInProbes;
  bool includeIndirectInProbes;
  Sample origSample;
  Renderer renderer;
  Scene scene;
  List<Point> surfacePoints;
  List<Spectrum> c_in;
}
