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

class TriangleFilter extends Filter {
  TriangleFilter(double xw, double yw)
      : super(xw, yw);

  double evaluate(double x, double y) {
    return max(0.0, xWidth - x.abs()) * max(0.0, yWidth - y.abs());
  }

  static TriangleFilter Create(ParamSet ps) {
    // Find common filter parameters
    double xw = ps.findOneFloat('xwidth', 2.0);
    double yw = ps.findOneFloat('ywidth', 2.0);
    return new TriangleFilter(xw, yw);
  }
}
