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
 * 4x4 matrix used for transformations.
 */
class Matrix4x4 {
  final Float32List m;

  /**
   * Static identity matrix.
   */
  static final Matrix4x4 IDENTITY = new Matrix4x4();

  /**
   * Creates a new identity matrix.
   */
  Matrix4x4() :
    m = new Float32List(16) {
    m[0] = 1.0;
    m[5] = 1.0;
    m[10] = 1.0;
    m[15] = 1.0;
  }

  /**
   * Create a matrix as a copy of [other].
   */
  Matrix4x4.from(Matrix4x4 other) :
    m = new Float32List.fromList(other.m);

  /**
   * Create a matrix as a copy of [other].
   */
  Matrix4x4.fromList(List other) :
    m = new Float32List.fromList(other);

  /**
   * Create a matrix from a set of 16 numbers.
   */
  Matrix4x4.values(double m11, double m12, double m13, double m14,
                   double m21, double m22, double m23, double m24,
                   double m31, double m32, double m33, double m34,
                   double m41, double m42, double m43, double m44) :
    m = new Float32List(16) {
    m[0] = m11;
    m[1] = m12;
    m[2] = m13;
    m[3] = m14;

    m[4] = m21;
    m[5] = m22;
    m[6] = m23;
    m[7] = m24;

    m[8] = m31;
    m[9] = m32;
    m[10] = m33;
    m[11] = m34;

    m[12] = m41;
    m[13] = m42;
    m[14] = m43;
    m[15] = m44;
  }

  bool operator==(Matrix4x4 m2) {
    for (int i = 0; i < 16; ++i) {
      if (m[i] != m2.m[i]) {
        return false;
      }
    }
    return true;
  }

  /**
   * Copy the values from another matrix into this matrix.
   */
  Matrix4x4 copy(Matrix4x4 other) {
    m[0] = other.m[0];
    m[1] = other.m[1];
    m[2] = other.m[2];
    m[3] = other.m[3];
    m[4] = other.m[4];
    m[5] = other.m[5];
    m[6] = other.m[6];
    m[7] = other.m[7];
    m[8] = other.m[8];
    m[9] = other.m[9];
    m[10] = other.m[10];
    m[11] = other.m[11];
    m[12] = other.m[12];
    m[13] = other.m[13];
    m[14] = other.m[14];
    m[15] = other.m[15];
    return this;
  }

  /**
   * Create a new matrix that is a clone of this matrix.
   */
  Matrix4x4 clone() => new Matrix4x4.from(this);

  /**
   * Set the elements of the matrix from a list.
   */
  Matrix4x4 setFromList(List other) {
    m[0] = other[0];
    m[1] = other[1];
    m[2] = other[2];
    m[3] = other[3];
    m[4] = other[4];
    m[5] = other[5];
    m[6] = other[6];
    m[7] = other[7];
    m[8] = other[8];
    m[9] = other[9];
    m[10] = other[10];
    m[11] = other[11];
    m[12] = other[12];
    m[13] = other[13];
    m[14] = other[14];
    m[15] = other[15];
    return this;
  }

  /**
   * Set the elements of the matrix.
   */
  Matrix4x4 set(double m11, double m12, double m13, double m14,
           double m21, double m22, double m23, double m24,
           double m31, double m32, double m33, double m34,
           double m41, double m42, double m43, double m44) {
    m[0] = m11;
    m[1] = m12;
    m[2] = m13;
    m[3] = m14;

    m[4] = m21;
    m[5] = m22;
    m[6] = m23;
    m[7] = m24;

    m[8] = m31;
    m[9] = m32;
    m[10] = m33;
    m[11] = m34;

    m[12] = m41;
    m[13] = m42;
    m[14] = m43;
    m[15] = m44;
    return this;
  }

  /**
   * Convert the matrix to a Float32List of 16 floats.
   */
  Float32List toFloat32List() => m;

  /**
   * Access a matrix element with the array operator.
   * Indexes are stored in row dominant order:
   * 0  1  2  3
   * 4  5  6  7
   * 8  9  10 11
   * 12 13 14 15
   */
  double operator [](int index) => m[index];

  /**
   * Set an element in the matrix.
   */
  operator []=(int index, double value) => m[index] = value;

  static Matrix4x4 Transpose(Matrix4x4 m) {
    return new Matrix4x4.values(m.m[0], m.m[4], m.m[8], m.m[12],
                                m.m[1], m.m[5], m.m[9], m.m[13],
                                m.m[2], m.m[6], m.m[10], m.m[14],
                                m.m[3], m.m[7], m.m[11], m.m[15]);
  }

  static Matrix4x4 Mul(Matrix4x4 m1, Matrix4x4 m2) {
    Matrix4x4 r = new Matrix4x4();

    for (int i = 0, k = 0; i < 4; ++i, k += 4) {
      for (int j = 0; j < 4; ++j) {
        r.m[k + j] = m1.m[k] * m2.m[j] +
                     m1.m[k + 1] * m2.m[4 + j] +
                     m1.m[k + 2] * m2.m[8 + j] +
                     m1.m[k + 3] * m2.m[12 + j];
      }
    }

    return r;
  }

  static Matrix4x4 Inverse(Matrix4x4 m) {
    return new Matrix4x4.from(m).invert();
  }

  /**
   * Set the matrix to the identity matrix.
   */
  Matrix4x4 setIdentity() {
    m[0] = 1.0;
    m[1] = 0.0;
    m[2] = 0.0;
    m[3] = 0.0;

    m[4] = 0.0;
    m[5] = 1.0;
    m[6] = 0.0;
    m[7] = 0.0;

    m[8] = 0.0;
    m[9] = 0.0;
    m[10] = 1.0;
    m[11] = 0.0;

    m[12] = 0.0;
    m[13] = 0.0;
    m[14] = 0.0;
    m[15] = 1.0;
    return this;
  }

  double determinant() {
    final double n11 = m[0];
    final double n12 = m[4];
    final double n13 = m[8];
    final double n14 = m[12];
    final double n21 = m[1];
    final double n22 = m[5];
    final double n23 = m[9];
    final double n24 = m[13];
    final double n31 = m[2];
    final double n32 = m[6];
    final double n33 = m[10];
    final double n34 = m[14];
    final double n41 = m[3];
    final double n42 = m[7];
    final double n43 = m[11];
    final double n44 = m[15];

    return
      (n14 * n23 * n32 * n41) -
      (n13 * n24 * n32 * n41) -
      (n14 * n22 * n33 * n41) +
      (n12 * n24 * n33 * n41) +

      (n13 * n22 * n34 * n41) -
      (n12 * n23 * n34 * n41) -
      (n14 * n23 * n31 * n42) +
      (n13 * n24 * n31 * n42) +

      (n14 * n21 * n33 * n42) -
      (n11 * n24 * n33 * n42) -
      (n13 * n21 * n34 * n42) +
      (n11 * n23 * n34 * n42) +

      (n14 * n22 * n31 * n43) -
      (n12 * n24 * n31 * n43) -
      (n14 * n21 * n32 * n43) +
      (n11 * n24 * n32 * n43) +

      (n12 * n21 * n34 * n43) -
      (n11 * n22 * n34 * n43) -
      (n13 * n22 * n31 * n44) +
      (n12 * n23 * n31 * n44) +

      (n13 * n21 * n32 * n44) -
      (n11 * n23 * n32 * n44) -
      (n12 * n21 * n33 * n44) +
      (n11 * n22 * n33 * n44);
  }

  /**
   * Invert the matrix
   */
  Matrix4x4 invert() {
    final double det = determinant();
    if (det == 0.0) {
      return this;
    }

    final double n11 = m[0];
    final double n12 = m[4];
    final double n13 = m[8];
    final double n14 = m[12];
    final double n21 = m[1];
    final double n22 = m[5];
    final double n23 = m[9];
    final double n24 = m[13];
    final double n31 = m[2];
    final double n32 = m[6];
    final double n33 = m[10];
    final double n34 = m[14];
    final double n41 = m[3];
    final double n42 = m[7];
    final double n43 = m[11];
    final double n44 = m[15];

    final double invDet = 1.0 / det;

    m[0] = (n23 * n34 * n42 - n24 * n33 * n42 + n24 * n32 * n43 -
               n22 * n34 * n43 - n23 * n32 * n44 + n22 * n33 * n44) * invDet;
    m[4] = (n14 * n33 * n42 - n13 * n34 * n42 - n14 * n32 * n43 +
           n12 * n34 * n43 + n13 * n32 * n44 - n12 * n33 * n44) * invDet;
    m[8] = (n13 * n24 * n42 - n14 * n23 * n42 + n14 * n22 * n43 -
           n12 * n24 * n43 - n13 * n22 * n44 + n12 * n23 * n44) * invDet;
    m[12] = (n14 * n23 * n32 - n13 * n24 * n32 - n14 * n22 * n33 +
            n12 * n24 * n33 + n13 * n22 * n34 - n12 * n23 * n34) * invDet;
    m[1] = (n24 * n33 * n41 - n23 * n34 * n41 - n24 * n31 * n43 +
           n21 * n34 * n43 + n23 * n31 * n44 - n21 * n33 * n44) * invDet;
    m[5] = (n13 * n34 * n41 - n14 * n33 * n41 + n14 * n31 * n43 -
           n11 * n34 * n43 - n13 * n31 * n44 + n11 * n33 * n44) * invDet;
    m[9] = (n14 * n23 * n41 - n13 * n24 * n41 - n14 * n21 * n43 +
           n11 * n24 * n43 + n13 * n21 * n44 - n11 * n23 * n44) * invDet;
    m[13] = (n13 * n24 * n31 - n14 * n23 * n31 + n14 * n21 * n33 -
            n11 * n24 * n33 - n13 * n21 * n34 + n11 * n23 * n34) * invDet;
    m[2] = (n22 * n34 * n41 - n24 * n32 * n41 + n24 * n31 * n42 -
           n21 * n34 * n42 - n22 * n31 * n44 + n21 * n32 * n44) * invDet;
    m[6] = (n14 * n32 * n41 - n12 * n34 * n41 - n14 * n31 * n42 +
           n11 * n34 * n42 + n12 * n31 * n44 - n11 * n32 * n44) * invDet;
    m[10] = (n12 * n24 * n41 - n14 * n22 * n41 + n14 * n21 * n42 -
            n11 * n24 * n42 - n12 * n21 * n44 + n11 * n22 * n44) * invDet;
    m[14] = (n14 * n22 * n31 - n12 * n24 * n31 - n14 * n21 * n32 +
            n11 * n24 * n32 + n12 * n21 * n34 - n11 * n22 * n34) * invDet;
    m[3] = (n23 * n32 * n41 - n22 * n33 * n41 - n23 * n31 * n42 +
           n21 * n33 * n42 + n22 * n31 * n43 - n21 * n32 * n43) * invDet;
    m[7] = (n12 * n33 * n41 - n13 * n32 * n41 + n13 * n31 * n42 -
           n11 * n33 * n42 - n12 * n31 * n43 + n11 * n32 * n43) * invDet;
    m[11] = (n13 * n22 * n41 - n12 * n23 * n41 - n13 * n21 * n42 +
            n11 * n23 * n42 + n12 * n21 * n43 - n11 * n22 * n43) * invDet;
    m[15] = (n12 * n23 * n31 - n13 * n22 * n31 + n13 * n21 * n32 -
            n11 * n23 * n32 - n12 * n21 * n33 + n11 * n22 * n33) * invDet;

    return this;
  }

  /**
   * Convert the matrix to a string.
   */
  String toString() {
    return "${m[0]} ${m[4]} ${m[8]} ${m[12]} "
           "${m[1]} ${m[5]} ${m[9]} ${m[13]} "
           "${m[2]} ${m[6]} ${m[10]} ${m[14]} "
           "${m[3]} ${m[7]} ${m[11]} ${m[15]}";
  }
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
