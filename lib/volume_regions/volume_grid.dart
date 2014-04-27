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

class VolumeGridDensity extends DensityRegion {
  VolumeGridDensity(Spectrum sa, Spectrum ss, double gg,
                    Spectrum emit, this.extent, Transform v2w,
                    this.nx, this.ny, this.nz, List<double> d) :
    super(sa, ss, gg, emit, v2w) {
    _density = new Float64List.fromList(d);
  }

  static VolumeGridDensity Create(Transform volume2world, ParamSet params) {
    // Initialize common volume region parameters
    Spectrum sigma_a = params.findOneSpectrum('sigma_a', new Spectrum(0.0));
    Spectrum sigma_s = params.findOneSpectrum('sigma_s', new Spectrum(0.0));
    double g = params.findOneFloat('g', 0.0);
    Spectrum Le = params.findOneSpectrum('Le', new Spectrum(0.0));
    Point p0 = params.findOnePoint('p0', new Point(0.0, 0.0, 0.0));
    Point p1 = params.findOnePoint('p1', new Point(1.0, 1.0, 1.0));
    List<double> data = params.findFloat('density');
    if (data == null) {
      LogError('No \'density\' values provided for volume grid?');
      return null;
    }

    int nx = params.findOneInt('nx', 1);
    int ny = params.findOneInt('ny', 1);
    int nz = params.findOneInt('nz', 1);

    if (data.length != nx * ny * nz) {
      LogError('VolumeGridDensity has ${data.length} density values but '
               'nx*ny*nz = ${nx * ny * nz}');
      return null;
    }

    return new VolumeGridDensity(sigma_a, sigma_s, g, Le, new BBox(p0, p1),
                                 volume2world, nx, ny, nz, data);
  }

  BBox worldBound() {
    return Transform.Inverse(worldToVolume).transformBBox(extent);
  }

  bool intersectP(Ray r, List<double> t0, List<double> t1) {
    Ray ray = worldToVolume.transformRay(r);
    return extent.intersectP(ray, t0, t1);
  }

  double density(Point Pobj) {
    if (!extent.inside(Pobj)) {
      return 0.0;
    }

    // Compute voxel coordinates and offsets for _Pobj_
    Vector vox = extent.offset(Pobj);
    vox.x = vox.x * nx - 0.5;
    vox.y = vox.y * ny - 0.5;
    vox.z = vox.z * nz - 0.5;

    int vx = vox.x.floor();
    int vy = vox.y.floor();
    int vz = vox.z.floor();

    double dx = vox.x - vx;
    double dy = vox.y - vy;
    double dz = vox.z - vz;

    // Trilinearly interpolate density values to compute local density
    double d00 = Lerp(dx, d(vx, vy, vz),     d(vx+1, vy, vz));
    double d10 = Lerp(dx, d(vx, vy+1, vz),   d(vx+1, vy+1, vz));
    double d01 = Lerp(dx, d(vx, vy, vz+1),   d(vx+1, vy, vz+1));
    double d11 = Lerp(dx, d(vx, vy+1, vz+1), d(vx+1, vy+1, vz+1));
    double d0 = Lerp(dy, d00, d10);
    double d1 = Lerp(dy, d01, d11);

    return Lerp(dz, d0, d1);
  }

  double d(int x, int y, int z) {
    x = x.clamp(0, nx - 1);
    y = y.clamp(0, ny - 1);
    z = z.clamp(0, nz - 1);
    return _density[z * nx * ny + y * nx + x];
  }

  List<double> _density;
  int nx, ny, nz;
  BBox extent;
}
