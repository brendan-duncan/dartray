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

abstract class DensityRegion extends VolumeRegion {
  DensityRegion(Spectrum sa, Spectrum ss, double gg,
                Spectrum emit, Transform volumeToWorld)
      : sig_a = new Spectrum.from(sa),
        sig_s = new Spectrum.from(ss),
        le = new Spectrum.from(emit),
        g = gg,
        worldToVolume = Transform.Inverse(volumeToWorld);

  double density(Point Pobj);

  Spectrum sigma_a(Point p, Vector v, double d) {
    return sig_a * density(worldToVolume.transformPoint(p));
  }

  Spectrum sigma_s(Point p, Vector v, double d) {
    return sig_s * density(worldToVolume.transformPoint(p));
  }

  Spectrum sigma_t(Point p, Vector v, double d) {
    return (sig_a + sig_s) * (density(worldToVolume.transformPoint(p)));
  }

  Spectrum Lve(Point p, Vector v, double d) {
    return le * density(worldToVolume.transformPoint(p));
  }

  double p(Point p, Vector w, Vector wp, double d) {
    return PhaseHG(w, wp, g);
  }

  Spectrum tau(Ray r, [double stepSize = 1.0, double u = 0.5]) {
    List<double> t0 = [0.0];
    List<double> t1 = [0.0];
    double length = r.direction.length();
    if (length == 0.0) {
      return new Spectrum(0.0);
    }

    Ray rn = new Ray(r.origin, r.direction / length,
                     r.minDistance * length,
                     r.maxDistance * length,
                     r.time);
    if (!intersectP(rn, t0, t1)) {
      return new Spectrum(0.0);
    }

    Spectrum tau = new Spectrum(0.0);
    t0[0] += u * stepSize;
    while (t0[0] < t1[0]) {
      tau += sigma_t(rn.pointAt(t0[0]), -rn.direction, r.time);
      t0[0] += stepSize;
    }

    return tau * stepSize;
  }


  Spectrum sig_a;
  Spectrum sig_s;
  Spectrum le;
  double g;
  Transform worldToVolume;
}
