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

/**
 * A 3-dimensional point in space.
 */
class Point extends Vector {
  Point([double x = 0.0, double y = 0.0, double z = 0.0]) :
    super(x, y, z);

  Point.from(Vector other) :
    super.from(other);

  Point operator*(double s) =>
      new Point(x * s, y * s, z * s);

  Point operator/(double s) =>
      new Point(x / s, y / s, z / s);

  Point operator+(Vector p) =>
      new Point(x + p.x, y + p.y, z + p.z);

  Point operator-(Vector p) =>
      new Point(x - p.x, y - p.y, z - p.z);

  static double Distance(Point p1, Point p2) =>
      (p2 - p1).length();

  static double DistanceSquared(Point p1, Point p2) =>
      (p2 - p1).lengthSquared();
}

