/****************************************************************************
 *  Copyright (C) 2014 by Brendan Duncan.                                   *
 *                                                                          *
 *  This file is part of DartRay.                                           *
 *                                                                          *
 *  Licensed under the Apache License, Version 2.0 (the "License");         *
 *  you may not use this file except in compliance with the License.        *
 *  You may obtain a copy of the License at                                 *
 *                                                                          *
 *  http://www.apache.org/licenses/LICENSE-2.0                              *
 *                                                                          *
 *  Unless required by applicable law or agreed to in writing, software     *
 *  distributed under the License is distributed on an "AS IS" BASIS,       *
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 *  See the License for the specific language governing permissions and     *
 *  limitations under the License.                                          *
 *                                                                          *
 *   This project is based on PBRT v2 ; see http://www.pbrt.org             *
 *   pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.*
 ****************************************************************************/
part of core;

abstract class DensityRegion extends VolumeRegion {
  DensityRegion(RGBColor sa, RGBColor ss, double gg,
                RGBColor emit, Transform volumeToWorld) :
    sig_a = new RGBColor.from(sa),
    sig_s = new RGBColor.from(ss),
    le = new RGBColor.from(emit),
    g = gg,
    worldToVolume = Transform.Inverse(volumeToWorld);

  double density(Point Pobj);

  RGBColor sigma_a(Point p, Vector v, double d) {
    return sig_a.scaled(density(worldToVolume.transformPoint(p)));
  }

  RGBColor sigma_s(Point p, Vector v, double d) {
    return sig_s.scaled(density(worldToVolume.transformPoint(p)));
  }

  RGBColor sigma_t(Point p, Vector v, double d) {
    return (sig_a + sig_s).scale(density(worldToVolume.transformPoint(p)));
  }

  RGBColor Lve(Point p, Vector v, double d) {
    return le.scaled(density(worldToVolume.transformPoint(p)));
  }

  double p(Point p, Vector w, Vector wp, double d) {
    return PhaseHG(w, wp, g);
  }

  RGBColor tau(Ray r, [double stepSize = 1.0, double u = 0.5]) {
    List<double> t0 = [0.0];
    List<double> t1 = [0.0];
    double length = r.direction.length();
    if (length == 0.0) {
      return new RGBColor(0.0);
    }

    Ray rn = new Ray(r.origin, r.direction / length,
                     r.minDistance * length,
                     r.maxDistance * length,
                     r.time);
    if (!intersectP(rn, t0, t1)) {
      return new RGBColor(0.0);
    }

    RGBColor tau = new RGBColor(0.0);
    t0[0] += u * stepSize;
    while (t0[0] < t1[0]) {
      tau += sigma_t(rn.pointAt(t0[0]), -rn.direction, r.time);
      t0[0] += stepSize;
    }

    return tau.scaled(stepSize);
  }


  RGBColor sig_a;
  RGBColor sig_s;
  RGBColor le;
  double g;
  Transform worldToVolume;
}
