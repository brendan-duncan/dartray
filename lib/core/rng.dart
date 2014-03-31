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
 * Generates random numbers using the Mersenne Twister algorithm.
 * This is used instead of Math.Random so that results can be verified against
 * the c++ implementation of PBRT.
 */
class RNG {
  static bool UseMathRandom = true;

  RNG([int seed = 5489]) {
    if (UseMathRandom) {
      random = new Math.Random(seed);
    } else {
      mti = N + 1; // mti==N+1 means mt[N] is not initialized
      this.seed(seed);
    }
  }

  void seed(int seed) {
    if (UseMathRandom) {
      random = new Math.Random(seed);
    } else {
      mt[0] = seed & 0xffffffff;
      for (mti = 1; mti < N; mti++) {
        mt[mti] =
          (1812433253 * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
          // See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
          // In the previous versions, MSBs of the seed affect
          // only MSBs of the array mt[].
          // 2002/01/09 modified by Makoto Matsumoto
        mt[mti] &= 0xffffffff; // for >32 bit machines
      }
    }
  }

  double randomFloat() {
    Stats.RNG_STARTED_RANDOM_FLOAT();
    if (UseMathRandom) {
      double v = random.nextDouble();
      Stats.RNG_FINISHED_RANDOM_FLOAT();
      return v;
    }
    double v = (randomUInt() & 0xffffff) / 16777216.0;
    Stats.RNG_FINISHED_RANDOM_FLOAT();
    return v;
  }

  int randomUInt() {
    if (UseMathRandom) {
      return random.nextInt(0xffffffff);
    }

    int y;
    const int M = 397;
    const int MATRIX_A = 0x9908b0df; // constant vector a
    const int UPPER_MASK = 0x80000000; // most significant w-r bits
    const int LOWER_MASK = 0x7fffffff; // least significant r bits

    const List<int> mag01 = const [0x0, MATRIX_A];
    // mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { // generate N words at one time
      Stats.RNG_STARTED_TABLEGEN();
      int kk;

      if (mti == N + 1) {   // if Seed() has not been called,
        seed(5489); // default initial seed
      }

      for (kk = 0; kk < N - M; kk++) {
        y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
        mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
      }

      for (; kk < N - 1;kk++) {
        y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
        mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
      }

      y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];

      mti = 0;
      Stats.RNG_FINISHED_TABLEGEN();
    }

    y = mt[mti++];

    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680;
    y ^= (y << 15) & 0xefc60000;
    y ^= (y >> 18);

    return y;
  }

  Math.Random random;

  static const int N = 624;
  Uint32List mt = new Uint32List(N); // the array for the state vector
  int mti;
}
