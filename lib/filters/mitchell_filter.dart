/****************************************************************************
 * Copyright (C) 2014 by Brendan Duncan.                                    *
 *                                                                          *
 * This file is part of DartRay.                                            *
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the 'License');          *
 * you may not use this file except in compliance with the License.         *
 * You may obtain a copy of the License at                                  *
 *                                                                          *
 * http://www.apache.org/licenses/LICENSE-2.0                               *
 *                                                                          *
 * Unless required by applicable law or agreed to in writing, software      *
 * distributed under the License is distributed on an 'AS IS' BASIS,        *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 * See the License for the specific language governing permissions and      *
 * limitations under the License.                                           *
 *                                                                          *
 * This project is based on PBRT v2 ; see http://www.pbrt.org               *
 * pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.  *
 ****************************************************************************/
part of filters;

class MitchellFilter extends Filter {
  MitchellFilter(this.b, this.c, double xw, double yw)
      : super(xw, yw);

  double evaluate(double x, double y) {
    return _mitchell1D(x * invXWidth) * _mitchell1D(y * invYWidth);
  }

  double _mitchell1D(double x) {
    x = (2.0 * x).abs();
    if (x > 1.0) {
      return ((-b - 6 * c) * x * x * x + (6 * b + 30 * c) * x * x +
              (-12 * b - 48 * c) * x + (8 * b + 24 * c)) * (1.0 / 6.0);
    } else {
      return ((12 - 9 * b - 6 * c) * x * x * x +
              (-18 + 12 * b + 6 * c) * x * x +
              (6 - 2 * b)) * (1.0 / 6.0);
    }
  }

  static MitchellFilter Create(ParamSet ps) {
    double xw = ps.findOneFloat('xwidth', 2.0);
    double yw = ps.findOneFloat('ywidth', 2.0);
    double B = ps.findOneFloat('B', 1.0 / 3.0);
    double C = ps.findOneFloat('C', 1.0 / 3.0);
    return new MitchellFilter(B, C, xw, yw);
  }

  double b;
  double c;
}
