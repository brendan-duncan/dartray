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

class AggregateVolume extends VolumeRegion {
  List<VolumeRegion> regions = [];

  AggregateVolume(this.regions)
      : _bound = new BBox() {
    for (int i = 0; i < regions.length; ++i) {
      _bound = BBox.Union(_bound, regions[i].worldBound());
    }
  }

  BBox worldBound() {
    return _bound;
  }

  bool intersectP(Ray ray, List<double> t0, List<double> t1) {
    t0[0] = INFINITY;
    t1[0] = -INFINITY;
    List<double> tr0 = [0.0];
    List<double> tr1 = [0.0];
    for (int i = 0; i < regions.length; ++i) {
      if (regions[i].intersectP(ray, tr0, tr1)) {
        t0[0] = Math.min(t0[0], tr0[0]);
        t1[0] = Math.max(t1[0], tr1[0]);
      }
    }
    return (t0[0] < t1[0]);
  }

  Spectrum sigma_a(Point p, Vector w, double time) {
    Spectrum s = new Spectrum(0.0);
    for (int i = 0; i < regions.length; ++i) {
      s += regions[i].sigma_a(p, w, time);
    }
    return s;
  }

  Spectrum sigma_s(Point p, Vector w, double time) {
    Spectrum s = new Spectrum(0.0);
    for (int i = 0; i < regions.length; ++i) {
      s += regions[i].sigma_s(p, w, time);
    }
    return s;
  }

  Spectrum Lve(Point p, Vector w, double time) {
    Spectrum L = new Spectrum(0.0);
    for (int i = 0; i < regions.length; ++i) {
      L += regions[i].Lve(p, w, time);
    }
    return L;
  }

  double p(Point p, Vector w, Vector wp, double time) {
    double ph = 0.0;
    double sumWt = 0.0;
    for (int i = 0; i < regions.length; ++i) {
      double wt = regions[i].sigma_s(p, w, time).luminance();
      sumWt += wt;
      ph += wt * regions[i].p(p, w, wp, time);
    }
    return ph / sumWt;
  }

  Spectrum sigma_t(Point p, Vector w, double time) {
    Spectrum s = new Spectrum(0.0);
    for (int i = 0; i < regions.length; ++i) {
      s += regions[i].sigma_t(p, w, time);
    }
    return s;
  }

  Spectrum tau(Ray ray, [double step = 1.0, double offset = 0.5]) {
    Spectrum t = new Spectrum(0.0);
    for (int i = 0; i < regions.length; ++i) {
      t += regions[i].tau(ray, step, offset);
    }
    return t;
  }

  BBox _bound;
}
