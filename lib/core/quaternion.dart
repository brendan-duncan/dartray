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
 * Represents a rotation in a suitable way for interpolation.
 */
class Quaternion {
  Vector v;
  double w;

  Quaternion()
      : v = new Vector(),
        w = 1.0;

  Quaternion.from(Quaternion q)
      : v = new Vector.from(q.v),
        w = q.w;

  Quaternion.fromMatrix(Matrix4x4 m)
      : v = new Vector(),
        w = 1.0 {
    double trace = m.data[0] + m.data[5] + m.data[10];
    if (trace > 0.0) {
      // Compute w from matrix trace, then xyz
      // 4w^2 = m[0][0] + m[1][1] + m[2][2] + m[3][3] (but m[3][3] == 1)
      double s = Math.sqrt(trace + 1.0);
      w = s / 2.0;
      s = 0.5 / s;
      v.x = (m.data[9] - m.data[6]) * s;
      v.y = (m.data[2] - m.data[8]) * s;
      v.z = (m.data[4] - m.data[1]) * s;
    } else {
      // Compute largest of $x$, $y$, or $z$, then remaining components
      const List<int> nxt = const [1, 2, 0];
      List<double> q = [0.0, 0.0, 0.0];
      int i = 0;
      if (m.data[5] > m.data[0]) {
        i = 1;
      }
      if (m.data[10] > m.data[i * 4 + i]) {
        i = 2;
      }
      int j = nxt[i];
      int k = nxt[j];
      double s = Math.sqrt((m.data[i * 4 + i] - (m.data[j * 4 + j] +
                           m.data[k * 4 + k])) + 1.0);
      q[i] = s * 0.5;
      if (s != 0.0) {
        s = 0.5 / s;
      }
      w = (m.data[k * 4 + j] - m.data[j * 4 + k]) * s;
      q[j] = (m.data[j * 4 + i] + m.data[i * 4 + j]) * s;
      q[k] = (m.data[k * 4 + i] + m.data[i * 4 + k]) * s;
      v.x = q[0];
      v.y = q[1];
      v.z = q[2];
    }
  }

  Quaternion copy(Quaternion other) {
    v = new Vector.from(other.v);
    w = other.w;
    return this;
  }

  Quaternion add(Quaternion q) {
    v.add(q.v);
    w += q.w;
    return this;
  }

  Quaternion operator +(Quaternion q) {
    return new Quaternion.from(this).add(q);
  }

  Quaternion sub(Quaternion q) {
    v.subtract(q.v);
    w -= q.w;
    return this;
  }

  Quaternion operator -(Quaternion q) {
    return new Quaternion.from(this).sub(q);
  }

  Quaternion scale(double f) {
    v.scale(f);
    w *= f;
    return this;
  }

  Quaternion operator *(double f) {
    return new Quaternion.from(this).scale(f);
  }

  Quaternion invScale(double f) {
    v.invScale(f);
    w /= f;
    return this;
  }

  Quaternion operator /(double f) {
    return new Quaternion.from(this).invScale(f);
  }

  Transform toTransform() {
    double xx = v.x * v.x;
    double yy = v.y * v.y;
    double zz = v.z * v.z;
    double xy = v.x * v.y;
    double xz = v.x * v.z;
    double yz = v.y * v.z;
    double wx = v.x * w;
    double wy = v.y * w;
    double wz = v.z * w;

    Matrix4x4 m = new Matrix4x4();
    m.data[0] = 1.0 - 2.0 * (yy + zz);
    m.data[1] =       2.0 * (xy + wz);
    m.data[2] =       2.0 * (xz - wy);
    m.data[4] =       2.0 * (xy - wz);
    m.data[5] = 1.0 - 2.0 * (xx + zz);
    m.data[6] =       2.0 * (yz + wx);
    m.data[8] =       2.0 * (xz + wy);
    m.data[9] =       2.0 * (yz - wx);
    m.data[10] = 1.0 - 2.0 * (xx + yy);

    // Transpose since we are left-handed.  Ugh.
    return new Transform(Matrix4x4.Transpose(m), m);
  }

  static Quaternion Slerp(double t, Quaternion q1, Quaternion q2) {
    double cosTheta = Dot(q1, q2);
    if (cosTheta > 0.9995) {
      return Normalize(q1 * (1.0 - t) + q2 * t);
    } else {
      double theta = Math.acos(cosTheta.clamp(-1.0, 1.0));
      double thetap = theta * t;
      Quaternion qperp = Normalize(q2 - (q1 * cosTheta));
      return q1 * Math.cos(thetap) + qperp * Math.sin(thetap);
    }
  }

  static double Dot(Quaternion q1, Quaternion q2) {
    return Vector.Dot(q1.v, q2.v) + q1.w * q2.w;
  }

  static Quaternion Normalize(Quaternion q) {
    return q / Math.sqrt(Dot(q, q));
  }
}
