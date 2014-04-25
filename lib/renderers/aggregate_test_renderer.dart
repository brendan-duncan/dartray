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

class AggregateTestRenderer extends Renderer {
  AggregateTestRenderer(this.nIterations, List<Primitive> p) {
    for (int i = 0; i < p.length; ++i) {
      p[i].fullyRefine(primitives);
    }

    for (int i = 0; i < primitives.length; ++i) {
      bboxes.add(primitives[i].worldBound());
    }
  }

  Future<OutputImage> render(Scene scene) {
    RNG rng = new RNG();
    // Compute bounding box of region used to generate random rays
    BBox bbox = new BBox.from(scene.worldBound);
    bbox.expand(bbox.pMax[bbox.maximumExtent()] -
                bbox.pMin[bbox.maximumExtent()]);

    Point lastHit;
    double lastEps = 0.0;
    for (int i = 0; i < nIterations; ++i) {
      // Choose random rays, _rayAccel_ and _rayAll_ for testing

      // Choose ray origin for testing accelerator
      Point org = new Point(Lerp(rng.randomFloat(), bbox.pMin.x, bbox.pMax.x),
                            Lerp(rng.randomFloat(), bbox.pMin.y, bbox.pMax.y),
                            Lerp(rng.randomFloat(), bbox.pMin.z, bbox.pMax.z));
      if ((rng.randomUint() % 4) == 0) {
        org = lastHit;
      }

      // Choose ray direction for testing accelerator
      Vector dir = UniformSampleSphere(rng.randomFloat(), rng.randomFloat());

      if ((rng.randomUint() % 32) == 0) {
        dir.x = dir.y = 0.0;
      } else if ((rng.randomUint() % 32) == 0) {
        dir.x = dir.z = 0.0;
      } else if ((rng.randomUint() % 32) == 0) {
        dir.y = dir.z = 0.0;
      }

      // Choose ray epsilon for testing accelerator
      double eps = 0.0;
      if (rng.randomFloat() < 0.25) {
        eps = lastEps;
      } else if (rng.randomFloat() < 0.25) {
        eps = 1.0e-3;
      }

      Ray rayAccel = new Ray(org, dir, eps);
      Ray rayAll = new Ray.from(rayAccel);

      // Compute intersections using accelerator and exhaustive testing
      Intersection isectAccel = new Intersection();
      Intersection isectAll = new Intersection();
      bool hitAccel = scene.intersect(rayAccel, isectAccel);
      bool hitAll = false;
      bool inconsistentBounds = false;

      for (int j = 0; j < primitives.length; ++j) {
        if (bboxes[j].intersectP(rayAll)) {
          hitAll = primitives[j].intersect(rayAll, isectAll) || hitAll;
        } else if (primitives[j].intersect(rayAll, isectAll)) {
          inconsistentBounds = true;
        }
      }

      // Report any inconsistencies between intersections
      if (!inconsistentBounds && ((hitAccel != hitAll) ||
          (rayAccel.maxDistance != rayAll.maxDistance))) {
        /*LogWarning("Disagreement: t accel %.16g [%a] t exhaustive %.16g [%a]\n"
                   "Ray: org [%a, %a, %a], dir [%a, %a, %a], mint = %a",
                   rayAccel.maxt, rayAll.maxt, rayAccel.maxt, rayAll.maxt,
                   rayAll.o.x, rayAll.o.y, rayAll.o.z,
                   rayAll.d.x, rayAll.d.y, rayAll.d.z, rayAll.mint);*/
      }

      if (hitAll) {
        lastHit = rayAll.pointAt(rayAll.maxDistance);
        lastEps = isectAll.rayEpsilon;
      }
    }

    Completer<OutputImage> c = new Completer<OutputImage>();
    c.complete(null);
    return c.future;
  }

  Spectrum Li(Scene scene, RayDifferential ray,
              Sample sample, RNG rng, [Intersection isect, Spectrum T]) {
    return new Spectrum(0.0);
  }

  Spectrum transmittance(Scene scene, RayDifferential ray,
          Sample sample, RNG rng) {
    return new Spectrum(0.0);
  }

  static AggregateTestRenderer Create(ParamSet params,
                                     List<Primitive> primitives) {
    int niters = params.findOneInt('niters', 100000);
    return new AggregateTestRenderer(niters, primitives);
  }

  int nIterations;
  List<Primitive> primitives = [];
  List<BBox> bboxes = [];
}
