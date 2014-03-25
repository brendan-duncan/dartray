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

const double INV_PI = 0.31830988618379067154;
const double INV_TWOPI = 0.15915494309189533577;
const double INV_FOURPI = 0.07957747154594766788;
const double INFINITY = 1.0e500;
const double FLT_EPSILON = 1.19209290e-07;

double Lerp(num t, num v1, num v2) =>
    (1.0 - t) * v1 + t * v2;

double Radians(num deg) =>
  (Math.PI / 180.0) * deg;

double Degrees(num rad) =>
  (180.0 / Math.PI) * rad;

double _invLog2 = 1.0 / Math.log(2.0);

double Log2(num x) {
  return Math.log(x) * _invLog2;
}

bool IsPowerOf2(int v) {
  return (v & (v - 1)) == 0;
}

int RoundUpPow2(int v) {
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  return v + 1;
}

int Mod(int a, int b) {
  int n = a ~/ b;
  a -= n * b;
  if (a < 0) {
    a += b;
  }
  return a;
}

double SmoothStep(double min, double max, double value) {
  double v = ((value - min) / (max - min)).clamp(0.0, 1.0);
  return v * v * (-2.0 * v  + 3.0);
}

bool Quadratic(double A, double B, double C, List<double> t0,
               List<double> t1) {
  // Find quadratic discriminant
  double discrim = B * B - 4.0 * A * C;
  if (discrim < 0.0) {
    return false;
  }

  double rootDiscrim = Math.sqrt(discrim);

  // Compute quadratic _t_ values
  double q;
  if (B < 0.0) {
    q = -0.5 * (B - rootDiscrim);
  } else {
    q = -0.5 * (B + rootDiscrim);
  }

  t0[0] = q / A;
  t1[0] = C / q;

  if (t0[0] > t1[0]) {
    double t = t0[0];
    t0[0] = t1[0];
    t1[0] = t;
  }

  return true;
}

/**
 * Standard compare function used for searches.
 */
bool less_than(a, b) => a.compareTo(b) < 0;

typedef bool BinaryPredicate(a, b);

typedef bool UnaryPredicate(a);

/**
 * Rearranges list such that all the elements which compare returns true
 * precede all those for which it returns false.  The index returned is
 * the starting element of the second group.
 */
int partition(List list, UnaryPredicate pred, [int first, int last]) {
  if (first == null) {
    first = 0;
  }
  if (last == null) {
    last = list.length;
  }

  while (first < last) {
    while (pred(list[first])) {
      ++first;
      if (first == last) {
        return first;
      }
    }

    do {
      --last;
      if (first == last) {
        return first;
      }
    } while (!pred(list[last]));

    var t = list[first];
    list[first] = list[last];
    list[last] = t;

    ++first;
  }
  return first;
}

void nth_element(List list, int first, int nth, int last,
                 BinaryPredicate pred) {
  // TODO implement nth_element properly for better performance.
  var l = list.sublist(first, last);
  l.sort((a, b) => pred(a, b) ? -1 : 1);
  for (int i = first, j = 0; i < last; ++i, ++j) {
    list[i] = l[j];
  }
}

/**
 * Returns an index of the first element in the sorted list which compares
 * greater than value.  Derived from c++ std::upper_bound.
 */
int upper_bound(List list, value, [BinaryPredicate compare = less_than]) {
  if (list.isEmpty) {
    return -1;
  }

  if (list.length == 1) {
    return 0;
  }

  int first = 0;
  int last = list.length - 1;
  int count = last;
  while (count > 0) {
    int index = first;
    int step = count >> 1;
    index += step;

    if (!compare(value, list[index])) {
      first = ++index;
      count -= step + 1;
    } else {
      count = step;
    }
  }

  return first;
}

/**
 * Returns an index of the first element in the sorted list which does not
 * compare less than value. Derived from c++ std::lower_bound.
 */
int lower_bound(List list, value, [BinaryPredicate compare = less_than]) {
  if (list.isEmpty) {
    return -1;
  }

  if (list.length == 1) {
    if (compare(list[0], value)) {
      return 1;
    }
    return 0;
  }

  int first = 0;
  int last = list.length - 1;
  int count = last;
  while (count > 0) {
    int index = first;
    int step = count ~/ 2;
    index += step;
    if (compare(list[index], value)) {
      first = ++index;
      count -= step + 1;
    } else {
      count = step;
    }
  }

  if (compare(list[first], value)) {
    return first + 1;
  }

  return first;
}

/**
 * Returns the index of the value in a sorted list, or -1 if the value is not
 * in the list. Derived from c++ std::binary_search.
 */
int find_element(List list, value, [bool compare(a, b) = less_than]) {
  if (list.isEmpty) {
    return -1;
  }

  if (list.length == 1) {
    if (compare(value, list[0]) == 0) {
      return 0;
    }
  }

  int index = lower_bound(list, value, compare);
  if (index != list.length) {
    if (!compare(value, list[index])) {
      return index;
    }
  }
  return -1;
}
