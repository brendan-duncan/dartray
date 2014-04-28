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

const double INV_PI = 0.31830988618379067154;
const double INV_TWOPI = 0.15915494309189533577;
const double INV_FOURPI = 0.07957747154594766788;
const double INFINITY = 1.0e500;
const double FLT_EPSILON = 1.19209290e-07;

/**
 * The Future class has a forEach static method, but not a while-loop
 * equivalent. This function will continue calling f while f returns a future,
 * and breaks when f returns null.
 * It waits for the future returned by f to complete before calling f again.
 * The returned future will complete when the loop has beeb broken.
 */
Future FutureWhileLoop(Function f) {
  Completer doneSignal = new Completer();
  void nextElement(_) {
    Future future = f();
    if (future == null) {
      doneSignal.complete(null);
    } else {
      new Future.sync(() => future)
        .then(nextElement, onError : doneSignal.completeError);
    }
  }
  nextElement(null);
  return doneSignal.future;
}

void GetSubWindow(int w, int h, int num, int count, List<int> extents) {
  // Determine how many tiles to use in each dimension, nx and ny
  int nx = count;
  int ny = 1;
  while ((nx & 0x1) == 0 && 2 * w * ny < h * nx) {
    nx >>= 1;
    ny <<= 1;
  }
  assert(nx * ny == count);

  // Compute x and y pixel sample range for sub-window
  int xo = num % nx;
  int yo = num ~/ nx;
  double tx0 = xo / nx;
  double tx1 = (xo + 1) / nx;
  double ty0 = yo / ny;
  double ty1 = (yo + 1) / ny;
  extents[0] = Lerp(tx0, 0, w).floor();
  extents[1] = Math.min(Lerp(tx1, 0, w).floor(), w);
  extents[2] = Lerp(ty0, 0, h).floor();
  extents[3] = Math.min(Lerp(ty1, 0, h).floor(), h);
}

/**
 * Linear interpolation between two values [v1] and [v2], at [t] which should
 * be between 0 and 1.
 */
Lerp(num t, v1, v2) =>
    v1 * (1.0 - t) + v2 * t;

/**
 * Convert degrees to radians.
 */
double Radians(num deg) =>
  (Math.PI / 180.0) * deg;

/**
 * Convert radians to degrees.
 */
double Degrees(num rad) =>
  (180.0 / Math.PI) * rad;

final double _invLog2 = 1.0 / Math.log(2.0);

double Log2(num x) {
  return Math.log(x) * _invLog2;
}

/**
 * Is the value [v] a power of 2?
 */
bool IsPowerOf2(int v) {
  return (v & (v - 1)) == 0;
}

/**
 * Round a value [v] to the next power of 2.
 */
int RoundUpPow2(int v) {
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  return v + 1;
}

/**
 * Smooth interpolation of [value] between [min] and [max].
 */
double SmoothStep(double min, double max, double value) {
  double v = ((value - min) / (max - min)).clamp(0.0, 1.0);
  return v * v * (-2.0 * v + 3.0);
}

/**
 * Solve the given quadratic equation.
 */
bool Quadratic(double A, double B, double C, List<double> t0, List<double> t1) {
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

bool SolveLinearSystem2x2(List<double> A, List<double> B,
                          List<double> x0, List<double> x1) {
  double det = A[0] * A[3] - A[1] * A[2];
  if (det.abs() < 1.0e-10) {
    return false;
  }

  x0[0] = (A[3] * B[0] - A[1] * B[1]) / det;
  x1[0] = (A[0] * B[1] - A[2] * B[0]) / det;

  if (x0[0].isNaN || x1[0].isNaN) {
    return false;
  }

  return true;
}

List<double> ReadFloatFile(List<int> bytes, String path) {
  String text = new String.fromCharCodes(bytes);
  int len = text.length;
  int ci = 0;

  final int ZERO = '0'.codeUnits[0];
  final int NINE = '9'.codeUnits[0];
  bool _isdigit(String c) {
    int cu = c.codeUnits[0];
    return cu >= ZERO && cu <= NINE;
  }

  bool _isspace(String c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r';
  }

  List<double> values = [];
  bool inNumber = false;
  String curNumber = '';
  int lineNumber = 0;
  while (ci < len) {
    String c = text[ci++];
    if (c == '\n') {
      ++lineNumber;
    }
    if (inNumber) {
      if (_isdigit(c) || c == '.' || c == 'e' || c == '-' || c == '+') {
        curNumber += c;
      } else {
        values.add(double.parse(curNumber));
        inNumber = false;
        curNumber = '';
      }
    } else {
      if (_isdigit(c) || c == '.' || c == '-' || c == '+') {
        inNumber = true;
        curNumber += c;
      } else if (c == '#') {
        while ((c = text[ci++]) != '\n' && ci < len);
        ++lineNumber;
      } else if (!_isspace(c)) {
        LogWarning('Unexpected text found at line $lineNumber of float file '
                   '$path: $c');
      }
    }
  }

  return values;
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
int upper_bound(List list, value, {int first: 0, int last,
                BinaryPredicate compare: less_than}) {
  if (list.isEmpty) {
    return -1;
  }

  if (list.length == 1) {
    return 0;
  }

  if (last == null) {
    last = list.length;
  }

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
int lower_bound(List list, value, {int first: 0, int last,
                BinaryPredicate compare: less_than}) {
  if (list.isEmpty) {
    return -1;
  }

  if (list.length == 1) {
    if (compare(list[0], value)) {
      return 1;
    }
    return 0;
  }

  if (last == null) {
    last = list.length;
  }

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
int find_element(List list, value, {int first: 0, int last,
                 BinaryPredicate compare: less_than}) {
  if (list.isEmpty) {
    return -1;
  }

  if (list.length == 1) {
    if (compare(value, list[0]) == 0) {
      return 0;
    }
  }

  if (last == null) {
    last = list.length;
  }

  int index = lower_bound(list, value, first: first, last: last,
                          compare: compare);

  if (index != list.length) {
    if (!compare(value, list[index])) {
      return index;
    }
  }
  return -1;
}

/**
 * This operation makes the elements in [first,last) into a heap.
 */
void make_heap(List list, int first, int last) {
  if (last - first < 2) {
    return;
  }

  int len = last - first;
  int parent = (len - 2) ~/ 2;
  while (true) {
    _adjust_heap(list, first, parent, len, list[first + parent]);
    if (parent == 0) {
      return;
    }
    parent--;
  }
}

void push_heap(List list, int first, int last) {
  _push_heap(list, first, (last - first) - 1, 0, list[last - 1]);
}

void pop_heap(List list, int first, int last) {
  _pop_heap(list, first, last - 1, last - 1, list[last - 1]);
}

void _adjust_heap(List list, int first, int holeIndex, int len, value) {
  final int topIndex = holeIndex;
  int secondChild = 2 * holeIndex + 2;
  while (secondChild < len) {
    if (list[first + secondChild] < list[first + (secondChild - 1)]) {
      secondChild--;
    }

    list[first + holeIndex] = list[first + secondChild];
    holeIndex = secondChild;
    secondChild = 2 * (secondChild + 1);
  }
  if (secondChild == len) {
    list[first + holeIndex] = list[first + (secondChild - 1)];
    holeIndex = secondChild - 1;
  }
  _push_heap(list, first, holeIndex, topIndex, value);
}

void _push_heap(List list, int first, int holeIndex, int topIndex, value) {
  int parent = (holeIndex - 1) ~/ 2;
  while (holeIndex > topIndex && list[first + parent] < value) {
    list[first + holeIndex] = list[first + parent];
    holeIndex = parent;
    parent = (holeIndex - 1) ~/ 2;
  }
  list[first + holeIndex] = value;
}

void _pop_heap(List list, int first, int last, int result, value) {
  list[result] = list[first];
  _adjust_heap(list, first, 0, last - first, value);
}
