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

class VisibilityTester {
  Ray r;

  void setSegment(Point p1, double eps1, Point p2, double eps2, double time) {
    double dist = Vector.Distance(p1, p2);
    r = new Ray(p1, (p2 - p1) / dist, eps1, dist * (1.0 - eps2), time);
  }

  void setRay(Point p, double eps, Vector w, double time) {
    r = new Ray(p, w, eps, INFINITY, time);
  }

  bool unoccluded(Scene scene) {
    return !scene.intersectP(r);
  }

  Spectrum transmittance(Scene scene, Renderer renderer, Sample sample,
                         RNG rng) {
    return renderer.transmittance(scene, new RayDifferential.fromRay(r),
                                  sample, rng);
  }
}
