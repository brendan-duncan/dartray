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

abstract class Texture {
  evaluate(DifferentialGeometry geom);
}

double Lanczos(double x, [double tau = 2.0]) {
  x = x.abs();
  if (x < 1.0e-5) {
    return 1.0;
  }
  if (x > 1.0) {
    return 0.0;
  }
  x *= Math.PI;
  double s = Math.sin(x * tau) / (x * tau);
  double lanczos = Math.sin(x) / x;
  return s * lanczos;
}

double Noise(double x, [double y = 0.5, double z = 0.5]) {
  // Compute noise cell coordinates and offsets
  int ix = x.floor();
  int iy = y.floor();
  int iz = z.floor();
  double dx = x - ix;
  double dy = y - iy;
  double dz = z - iz;

  // Compute gradient weights
  ix &= (_NOISE_PERM_SIZE - 1);
  iy &= (_NOISE_PERM_SIZE - 1);
  iz &= (_NOISE_PERM_SIZE - 1);

  double w000 = _Grad(ix,   iy,   iz,   dx,   dy,   dz);
  double w100 = _Grad(ix+1, iy,   iz,   dx-1, dy,   dz);
  double w010 = _Grad(ix,   iy+1, iz,   dx,   dy-1, dz);
  double w110 = _Grad(ix+1, iy+1, iz,   dx-1, dy-1, dz);
  double w001 = _Grad(ix,   iy,   iz+1, dx,   dy,   dz-1);
  double w101 = _Grad(ix+1, iy,   iz+1, dx-1, dy,   dz-1);
  double w011 = _Grad(ix,   iy+1, iz+1, dx,   dy-1, dz-1);
  double w111 = _Grad(ix+1, iy+1, iz+1, dx-1, dy-1, dz-1);

  // Compute trilinear interpolation of weights
  double wx = _NoiseWeight(dx);
  double wy = _NoiseWeight(dy);
  double wz = _NoiseWeight(dz);

  double x00 = Lerp(wx, w000, w100);
  double x10 = Lerp(wx, w010, w110);
  double x01 = Lerp(wx, w001, w101);
  double x11 = Lerp(wx, w011, w111);
  double y0 = Lerp(wy, x00, x10);
  double y1 = Lerp(wy, x01, x11);

  return Lerp(wz, y0, y1);
}

double NoisePoint(Point P) => Noise(P.x, P.y, P.z);

double FBm(Point P, Vector dpdx, Vector dpdy, double omega, int maxOctaves) {
  // Compute number of octaves for antialiased FBm
  double s2 = Math.max(dpdx.lengthSquared(), dpdy.lengthSquared());
  double log2_s2 = Log2(s2);
  double foctaves = Math.min(maxOctaves.toDouble(),
                             Math.max(0.0, -1.0 - 0.5 * log2_s2));
  int octaves = foctaves.floor();

  // Compute sum of octaves of noise for FBm
  double sum = 0.0;
  double lambda = 1.0;
  double o = 1.0;
  for (int i = 0; i < octaves; ++i) {
    sum += o * NoisePoint(P * lambda);
    lambda *= 1.99;
    o *= omega;
  }

  double partialOctave = foctaves - octaves;
  sum += o * SmoothStep(0.3, 0.7, partialOctave) * NoisePoint(P * lambda);

  return sum;
}

double Turbulence(Point P, Vector dpdx, Vector dpdy,
                  double omega, int maxOctaves) {
  // Compute number of octaves for antialiased FBm
  double s2 = Math.max(dpdx.lengthSquared(), dpdy.lengthSquared());
  double foctaves = Math.min(maxOctaves.toDouble(),
                             Math.max(0.0, -1.0 - 0.5 * Log2(s2)));
  int octaves = foctaves.floor();

  // Compute sum of octaves of noise for turbulence
  double sum = 0.0;
  double lambda = 1.0;
  double o = 1.0;
  for (int i = 0; i < octaves; ++i) {
    sum += o * NoisePoint(P * lambda).abs();
    lambda *= 1.99;
    o *= omega;
  }
  double partialOctave = foctaves - octaves;
  sum += o * SmoothStep(0.3, 0.7, partialOctave) *
         NoisePoint(P * lambda).abs();

  // finally, add in value to account for average value of fabsf(Noise())
  // (~0.2) for the remaining octaves...
  sum += (maxOctaves - foctaves) * 0.2;

  return sum;
}


const int _NOISE_PERM_SIZE = 256;

double _Grad(int x, int y, int z, double dx, double dy, double dz) {
  int h = _NOISE_PERM[_NOISE_PERM[_NOISE_PERM[x] + y] + z];
  h &= 15;
  double u = h < 8 || h == 12 || h == 13 ? dx : dy;
  double v = h < 4 || h == 12 || h == 13 ? dy : dz;
  return ((h & 1) != 0 ? -u : u) + ((h & 2) != 0 ? -v : v);
}


double _NoiseWeight(double t) {
  double t3 = t * t * t;
  double t4 = t3 * t;
  return 6.0 * t4 * t - 15.0 * t4 + 10.0 * t3;
}

const List<int> _NOISE_PERM = const[
  151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96,
  53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142,
  // Remainder of the noise permutation table
  8, 99, 37, 240, 21, 10, 23,
  190,  6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35,
  11, 32, 57, 177, 33,
  88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168,  68, 175, 74, 165, 71, 134,
  139, 48, 27, 166,
  77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55,
  46, 245, 40, 244,
  102, 143, 54,  65, 25, 63, 161,  1, 216, 80, 73, 209, 76, 132, 187, 208,  89,
  18, 169, 200, 196,
  135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186,  3, 64, 52, 217,
  226, 250, 124, 123,
  5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58,
  17, 182, 189, 28, 42,
  223, 183, 170, 213, 119, 248, 152,  2, 44, 154, 163,  70, 221, 153, 101, 155,
  167,  43, 172, 9,
  129, 22, 39, 253,  19, 98, 108, 110, 79, 113, 224, 232, 178, 185,  112, 104,
  218, 246, 97, 228,
  251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241,  81, 51, 145, 235,
  249, 14, 239, 107,
  49, 192, 214,  31, 181, 199, 106, 157, 184,  84, 204, 176, 115, 121, 50, 45,
  127,  4, 150, 254,
  138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215,
  61, 156, 180,
  151, 160, 137, 91, 90, 15,
  131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99,
  37, 240, 21, 10, 23,
  190,  6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35,
  11, 32, 57, 177, 33,
  88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168,  68, 175, 74, 165, 71, 134,
  139, 48, 27, 166,
  77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55,
  46, 245, 40, 244,
  102, 143, 54,  65, 25, 63, 161,  1, 216, 80, 73, 209, 76, 132, 187, 208,  89,
  18, 169, 200, 196,
  135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186,  3, 64, 52, 217,
  226, 250, 124, 123,
  5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58,
  17, 182, 189, 28, 42,
  223, 183, 170, 213, 119, 248, 152,  2, 44, 154, 163,  70, 221, 153, 101, 155,
  167,  43, 172, 9,
  129, 22, 39, 253,  19, 98, 108, 110, 79, 113, 224, 232, 178, 185,  112, 104,
  218, 246, 97, 228,
  251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241,  81, 51, 145, 235,
  249, 14, 239, 107,
  49, 192, 214,  31, 181, 199, 106, 157, 184,  84, 204, 176, 115, 121, 50, 45,
  127,  4, 150, 254,
  138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215,
  61, 156, 180
];
