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

/**
 * Animates objects and cameras over time as a linear interpolation between
 * two transforms.
 *
 * The transformation used to evaluate an object is linearly interpolated
 * between [transform1] at [startTime], and [transform2] at [endTime], for the
 * current camera shutter time being evaluated.
 *
 * If [transform1] and [transform2] are the same, then the transform is
 * considered not-animated, and the interpolation is short-cutted to return
 * [transform1].
 */
class AnimatedTransform {
  AnimatedTransform(Transform transform1, this.startTime,
                    Transform transform2, this.endTime)
      : startTransform = new Transform.from(transform1),
        endTransform = new Transform.from(transform2),
        actuallyAnimated = (transform1 != transform2) {
    Decompose(startTransform.m, T[0], R[0], S[0]);
    Decompose(endTransform.m, T[1], R[1], S[1]);
  }

  AnimatedTransform.from(AnimatedTransform other)
      : startTime = other.startTime,
        endTime = other.endTime,
        startTransform = new Transform.from(other.startTransform),
        endTransform = new Transform.from(other.endTransform),
        actuallyAnimated = other.actuallyAnimated {
    T[0] = new Vector.from(other.T[0]);
    T[1] = new Vector.from(other.T[1]);

    R[0] = new Quaternion.from(other.R[0]);
    R[1] = new Quaternion.from(other.R[1]);

    S[0] = new Matrix4x4.from(other.S[0]);
    S[1] = new Matrix4x4.from(other.S[1]);
  }

  static void Decompose(Matrix4x4 m, Vector T, Quaternion Rquat, Matrix4x4 S) {
    // Extract translation _T_ from transformation matrix
    T.x = m.data[3];
    T.y = m.data[7];
    T.z = m.data[11];

    // Compute new transformation matrix _M_ without translation
    Matrix4x4 M = new Matrix4x4.from(m);
    for (int i = 0; i < 3; ++i) {
      M.data[i * 4 + 3] = 0.0;
      M.data[12 + i] = 0.0;
    }
    M.data[15] = 1.0;

    // Extract rotation _R_ from transformation matrix
    double norm;
    int count = 0;
    Matrix4x4 R = new Matrix4x4.from(M);

    do {
      // Compute next matrix _Rnext_ in series
      Matrix4x4 Rnext = new Matrix4x4();
      Matrix4x4 Rit = Matrix4x4.Inverse(Matrix4x4.Transpose(R));
      for (int i = 0; i < 16; ++i) {
        Rnext.data[i] = 0.5 * (R.data[i] + Rit.data[i]);
      }

      // Compute norm of difference between _R_ and _Rnext_
      norm = 0.0;
      for (int i = 0, j = 0; i < 3; ++i, j += 4) {
        double n = (R.data[j] - Rnext.data[j]).abs() +
                   (R.data[j + 1] - Rnext.data[j + 1]).abs() +
                   (R.data[j + 2] - Rnext.data[j + 2]).abs();
        norm = Math.max(norm, n);
      }

      R = Rnext;
    } while (++count < 100 && norm > 0.0001);

    // XXX TODO FIXME deal with flip...
    Rquat.copy(new Quaternion.fromMatrix(R));

    // Compute scale S using rotation and original matrix
    S.copy(Matrix4x4.Mul(Matrix4x4.Inverse(R), M));
  }

  void interpolate(double time, Transform t) {
    // Handle boundary conditions for matrix interpolation
    if (!actuallyAnimated || time <= startTime) {
      t.copy(startTransform);
      return;
    }

    if (time >= endTime) {
      t.copy(endTransform);
      return;
    }

    double dt = (time - startTime) / (endTime - startTime);

    // Interpolate translation at dt
    Vector trans = T[0] * (1.0 - dt) + T[1] * dt;

    // Interpolate rotation at dt
    Quaternion rotate = Quaternion.Slerp(dt, R[0], R[1]);

    // Interpolate scale at dt
    Matrix4x4 scale = new Matrix4x4();
    for (int i = 0; i < 16; ++i) {
      scale.data[i] = Lerp(dt, S[0].data[i], S[1].data[i]);
    }

    // Compute interpolated matrix as product of interpolated components
    t.copy(Transform.Translate(trans) * rotate.toTransform() *
           new Transform(scale));
  }

  Ray transformRay(Ray r, [Ray tr]) {
    if (tr == null) {
      tr = new Ray();
    }

    if (!actuallyAnimated || r.time <= startTime) {
      startTransform.transformRay(r, tr);
    } else if (r.time >= endTime) {
      endTransform.transformRay(r, tr);
    } else {
      Transform t = new Transform();
      interpolate(r.time, t);
      t.transformRay(r, tr);
    }

    tr.time = r.time;

    return tr;
  }

  void transformRayDifferential(RayDifferential r, RayDifferential tr) {
    if (!actuallyAnimated || r.time <= startTime) {
      startTransform.transformRayDifferential(r, tr);
    } else if (r.time >= endTime) {
      endTransform.transformRayDifferential(r, tr);
    } else {
      Transform t = new Transform();
      interpolate(r.time, t);
      t.transformRayDifferential(r, tr);
    }
    tr.time = r.time;
  }

  Point transformPoint(double time, Point p) {
    if (!actuallyAnimated || time <= startTime) {
      return startTransform.transformPoint(p);
    } else if (time >= endTime) {
      return endTransform.transformPoint(p);
    }
    Transform t = new Transform();
    interpolate(time, t);
    return t.transformPoint(p);
  }

  Vector transformVector(double time, Vector v) {
    if (!actuallyAnimated || time <= startTime) {
      return startTransform.transformVector(v);
    } else if (time >= endTime) {
      return endTransform.transformVector(v);
    }
    Transform t = new Transform();
    interpolate(time, t);
    return t.transformVector(v);
  }

  BBox motionBounds(BBox b, bool useInverse) {
    if (!actuallyAnimated) {
      return Transform.Inverse(startTransform).transformBBox(b);
    }

    BBox ret = new BBox();
    const int nSteps = 128;
    Transform t = new Transform();
    for (int i = 0; i < nSteps; ++i) {
      double time = Lerp(i / (nSteps - 1), startTime, endTime);
      interpolate(time, t);
      if (useInverse) {
        t = Transform.Inverse(t);
      }
      ret = BBox.Union(ret, t.transformBBox(b));
    }
    return ret;
  }

  bool hasScale() {
    return startTransform.hasScale() || endTransform.hasScale();
  }

  double startTime;
  double endTime;
  Transform startTransform;
  Transform endTransform;
  bool actuallyAnimated;
  List<Vector> T = [new Vector(), new Vector()];
  List<Quaternion> R = [new Quaternion(), new Quaternion()];
  List<Matrix4x4> S = [new Matrix4x4(), new Matrix4x4()];
}
