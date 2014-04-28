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
 * 4x4 matrix used for transformations.
 */
class Matrix4x4 {
  final Float32List data;

  /**
   * Static identity matrix.
   */
  static final Matrix4x4 IDENTITY = new Matrix4x4();

  /**
   * Creates a new identity matrix.
   */
  Matrix4x4()
      : data = new Float32List(16) {
    data[0] = 1.0;
    data[5] = 1.0;
    data[10] = 1.0;
    data[15] = 1.0;
  }

  /**
   * Create a matrix as a copy of [other].
   */
  Matrix4x4.from(Matrix4x4 other)
      : data = new Float32List.fromList(other.data);

  /**
   * Create a matrix as a copy of [other].
   */
  Matrix4x4.fromList(List other)
      : data = new Float32List.fromList(other);

  /**
   * Create a matrix from a set of 16 numbers.
   */
  Matrix4x4.values(double m11, double m12, double m13, double m14,
                   double m21, double m22, double m23, double m24,
                   double m31, double m32, double m33, double m34,
                   double m41, double m42, double m43, double m44)
      : data = new Float32List(16) {
    data[0] = m11;
    data[1] = m12;
    data[2] = m13;
    data[3] = m14;

    data[4] = m21;
    data[5] = m22;
    data[6] = m23;
    data[7] = m24;

    data[8] = m31;
    data[9] = m32;
    data[10] = m33;
    data[11] = m34;

    data[12] = m41;
    data[13] = m42;
    data[14] = m43;
    data[15] = m44;
  }

  bool operator ==(Matrix4x4 m2) {
    for (int i = 0; i < 16; ++i) {
      if (data[i] != m2.data[i]) {
        return false;
      }
    }
    return true;
  }

  /**
   * Copy the values from another matrix into this matrix.
   */
  Matrix4x4 copy(Matrix4x4 other) {
    data[0] = other.data[0];
    data[1] = other.data[1];
    data[2] = other.data[2];
    data[3] = other.data[3];
    data[4] = other.data[4];
    data[5] = other.data[5];
    data[6] = other.data[6];
    data[7] = other.data[7];
    data[8] = other.data[8];
    data[9] = other.data[9];
    data[10] = other.data[10];
    data[11] = other.data[11];
    data[12] = other.data[12];
    data[13] = other.data[13];
    data[14] = other.data[14];
    data[15] = other.data[15];
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
    data[0] = other[0];
    data[1] = other[1];
    data[2] = other[2];
    data[3] = other[3];
    data[4] = other[4];
    data[5] = other[5];
    data[6] = other[6];
    data[7] = other[7];
    data[8] = other[8];
    data[9] = other[9];
    data[10] = other[10];
    data[11] = other[11];
    data[12] = other[12];
    data[13] = other[13];
    data[14] = other[14];
    data[15] = other[15];
    return this;
  }

  /**
   * Set the elements of the matrix.
   */
  Matrix4x4 set(double m11, double m12, double m13, double m14,
                double m21, double m22, double m23, double m24,
                double m31, double m32, double m33, double m34,
                double m41, double m42, double m43, double m44) {
    data[0] = m11;
    data[1] = m12;
    data[2] = m13;
    data[3] = m14;

    data[4] = m21;
    data[5] = m22;
    data[6] = m23;
    data[7] = m24;

    data[8] = m31;
    data[9] = m32;
    data[10] = m33;
    data[11] = m34;

    data[12] = m41;
    data[13] = m42;
    data[14] = m43;
    data[15] = m44;
    return this;
  }

  /**
   * Access a matrix element with the array operator.
   * Indexes are stored in row dominant order:
   * 0  1  2  3
   * 4  5  6  7
   * 8  9  10 11
   * 12 13 14 15
   */
  double operator [](int index) => data[index];

  /**
   * Set an element in the matrix.
   */
  operator []=(int index, double value) => data[index] = value;

  static Matrix4x4 Transpose(Matrix4x4 m) {
    return new Matrix4x4.values(m.data[0], m.data[4], m.data[8], m.data[12],
                                m.data[1], m.data[5], m.data[9], m.data[13],
                                m.data[2], m.data[6], m.data[10], m.data[14],
                                m.data[3], m.data[7], m.data[11], m.data[15]);
  }

  static Matrix4x4 Mul(Matrix4x4 m1, Matrix4x4 m2) {
    Matrix4x4 r = new Matrix4x4();

    for (int i = 0, k = 0; i < 4; ++i, k += 4) {
      for (int j = 0; j < 4; ++j) {
        r.data[k + j] = m1.data[k] * m2.data[j] +
                     m1.data[k + 1] * m2.data[4 + j] +
                     m1.data[k + 2] * m2.data[8 + j] +
                     m1.data[k + 3] * m2.data[12 + j];
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
    data[0] = 1.0;
    data[1] = 0.0;
    data[2] = 0.0;
    data[3] = 0.0;

    data[4] = 0.0;
    data[5] = 1.0;
    data[6] = 0.0;
    data[7] = 0.0;

    data[8] = 0.0;
    data[9] = 0.0;
    data[10] = 1.0;
    data[11] = 0.0;

    data[12] = 0.0;
    data[13] = 0.0;
    data[14] = 0.0;
    data[15] = 1.0;
    return this;
  }

  double determinant() {
    final double n11 = data[0];
    final double n12 = data[4];
    final double n13 = data[8];
    final double n14 = data[12];
    final double n21 = data[1];
    final double n22 = data[5];
    final double n23 = data[9];
    final double n24 = data[13];
    final double n31 = data[2];
    final double n32 = data[6];
    final double n33 = data[10];
    final double n34 = data[14];
    final double n41 = data[3];
    final double n42 = data[7];
    final double n43 = data[11];
    final double n44 = data[15];

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

    final double n11 = data[0];
    final double n12 = data[4];
    final double n13 = data[8];
    final double n14 = data[12];
    final double n21 = data[1];
    final double n22 = data[5];
    final double n23 = data[9];
    final double n24 = data[13];
    final double n31 = data[2];
    final double n32 = data[6];
    final double n33 = data[10];
    final double n34 = data[14];
    final double n41 = data[3];
    final double n42 = data[7];
    final double n43 = data[11];
    final double n44 = data[15];

    final double invDet = 1.0 / det;

    data[0] = (n23 * n34 * n42 - n24 * n33 * n42 + n24 * n32 * n43 -
               n22 * n34 * n43 - n23 * n32 * n44 + n22 * n33 * n44) * invDet;
    data[4] = (n14 * n33 * n42 - n13 * n34 * n42 - n14 * n32 * n43 +
               n12 * n34 * n43 + n13 * n32 * n44 - n12 * n33 * n44) * invDet;
    data[8] = (n13 * n24 * n42 - n14 * n23 * n42 + n14 * n22 * n43 -
               n12 * n24 * n43 - n13 * n22 * n44 + n12 * n23 * n44) * invDet;
    data[12] = (n14 * n23 * n32 - n13 * n24 * n32 - n14 * n22 * n33 +
                n12 * n24 * n33 + n13 * n22 * n34 - n12 * n23 * n34) * invDet;
    data[1] = (n24 * n33 * n41 - n23 * n34 * n41 - n24 * n31 * n43 +
               n21 * n34 * n43 + n23 * n31 * n44 - n21 * n33 * n44) * invDet;
    data[5] = (n13 * n34 * n41 - n14 * n33 * n41 + n14 * n31 * n43 -
               n11 * n34 * n43 - n13 * n31 * n44 + n11 * n33 * n44) * invDet;
    data[9] = (n14 * n23 * n41 - n13 * n24 * n41 - n14 * n21 * n43 +
               n11 * n24 * n43 + n13 * n21 * n44 - n11 * n23 * n44) * invDet;
    data[13] = (n13 * n24 * n31 - n14 * n23 * n31 + n14 * n21 * n33 -
                n11 * n24 * n33 - n13 * n21 * n34 + n11 * n23 * n34) * invDet;
    data[2] = (n22 * n34 * n41 - n24 * n32 * n41 + n24 * n31 * n42 -
               n21 * n34 * n42 - n22 * n31 * n44 + n21 * n32 * n44) * invDet;
    data[6] = (n14 * n32 * n41 - n12 * n34 * n41 - n14 * n31 * n42 +
               n11 * n34 * n42 + n12 * n31 * n44 - n11 * n32 * n44) * invDet;
    data[10] = (n12 * n24 * n41 - n14 * n22 * n41 + n14 * n21 * n42 -
                n11 * n24 * n42 - n12 * n21 * n44 + n11 * n22 * n44) * invDet;
    data[14] = (n14 * n22 * n31 - n12 * n24 * n31 - n14 * n21 * n32 +
                n11 * n24 * n32 + n12 * n21 * n34 - n11 * n22 * n34) * invDet;
    data[3] = (n23 * n32 * n41 - n22 * n33 * n41 - n23 * n31 * n42 +
               n21 * n33 * n42 + n22 * n31 * n43 - n21 * n32 * n43) * invDet;
    data[7] = (n12 * n33 * n41 - n13 * n32 * n41 + n13 * n31 * n42 -
               n11 * n33 * n42 - n12 * n31 * n43 + n11 * n32 * n43) * invDet;
    data[11] = (n13 * n22 * n41 - n12 * n23 * n41 - n13 * n21 * n42 +
                n11 * n23 * n42 + n12 * n21 * n43 - n11 * n22 * n43) * invDet;
    data[15] = (n12 * n23 * n31 - n13 * n22 * n31 + n13 * n21 * n32 -
                n11 * n23 * n32 - n12 * n21 * n33 + n11 * n22 * n33) * invDet;

    return this;
  }

  /**
   * Convert the matrix to a string.
   */
  String toString() {
    return '${data[0]} ${data[1]} ${data[2]} ${data[3]} '
           '${data[4]} ${data[5]} ${data[6]} ${data[7]} '
           '${data[8]} ${data[9]} ${data[10]} ${data[11]} '
           '${data[12]} ${data[13]} ${data[14]} ${data[15]}';
  }
}

