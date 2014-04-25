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
part of volume_regions;

class HomogeneousVolumeRegion extends VolumeRegion {
  HomogeneousVolumeRegion(Spectrum sa, Spectrum ss, double gg,
                                 Spectrum emit, BBox e, Transform v2w) {
    worldToVolume = Transform.Inverse(v2w);
    sig_a = sa;
    sig_s = ss;
    g = gg;
    le = emit;
    extent = e;
  }

  static HomogeneousVolumeRegion Create(Transform volume2world,
                                         ParamSet params) {
    // Initialize common volume region parameters
    Spectrum sigma_a = params.findOneSpectrum("sigma_a", new Spectrum(0.0));
    Spectrum sigma_s = params.findOneSpectrum("sigma_s", new Spectrum(0.0));
    double g = params.findOneFloat("g", 0.0);
    Spectrum Le = params.findOneSpectrum("Le", new Spectrum(0.0));
    Point p0 = params.findOnePoint("p0", new Point(0.0, 0.0, 0.0));
    Point p1 = params.findOnePoint("p1", new Point(1.0, 1.0, 1.0));

    return new HomogeneousVolumeRegion(sigma_a, sigma_s, g, Le,
                                        new BBox(p0, p1),
                                        volume2world);
  }

  BBox worldBound() {
    return Transform.Inverse(worldToVolume).transformBBox(extent);
  }

  bool intersectP(Ray r, List<double> t0, List<double> t1) {
    Ray ray = worldToVolume.transformRay(r);
    return extent.intersectP(ray, t0, t1);
  }

  Spectrum sigma_a(Point p, Vector v, double d) {
    return extent.inside(worldToVolume.transformPoint(p)) ? sig_a :
           new Spectrum(0.0);
  }

  Spectrum sigma_s(Point p, Vector v, double d) {
      return extent.inside(worldToVolume.transformPoint(p)) ? sig_s :
             new Spectrum(0.0);
  }

  Spectrum sigma_t(Point p, Vector v, double d) {
      return extent.inside(worldToVolume.transformPoint(p)) ? (sig_a + sig_s) :
             new Spectrum(0.0);
  }

  Spectrum Lve(Point p, Vector v, double d) {
      return extent.inside(worldToVolume.transformPoint(p)) ? le :
             new Spectrum(0.0);
  }

  double p(Point p, Vector wi, Vector wo, double time) {
    if (!extent.inside(worldToVolume.transformPoint(p))) {
      return 0.0;
    }
    return PhaseHG(wi, wo, g);
  }

  Spectrum tau(Ray ray, [double step = 1.0, double offset = 0.5]) {
    List<double> t0 = [0.0];
    List<double> t1 = [0.0];
    if (!intersectP(ray, t0, t1)) {
      return new Spectrum(0.0);
    }
    return (sig_a + sig_s) *
           Vector.Distance(ray.pointAt(t0[0]), ray.pointAt(t1[0]));
  }

  Spectrum sig_a, sig_s, le;
  double g;
  BBox extent;
  Transform worldToVolume;
}
