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

class LanczosSincFilter extends Filter {
  LanczosSincFilter(double xw, double yw, this.tau)
      : super(xw, yw);

  double evaluate(double x, double y) {
    return _sinc1D(x * invXWidth) * _sinc1D(y * invYWidth);
  }

  double _sinc1D(double x) {
    x = x.abs();
    if (x < 1e-5) {
      return 1.0;
    }
    if (x > 1.0) {
      return 0.0;
    }
    x *= PI;
    double sinc = sin(x) / x;
    double lanczos = sin(x * tau) / (x * tau);
    return sinc * lanczos;
  }

  static LanczosSincFilter Create(ParamSet ps) {
    double xw = ps.findOneFloat('xwidth', 4.0);
    double yw = ps.findOneFloat('ywidth', 4.0);
    double tau = ps.findOneFloat('tau', 3.0);
    return new LanczosSincFilter(xw, yw, tau);
  }

  double tau;
}
