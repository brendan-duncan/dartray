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
part of pixel_samplers;

class TilePixelSampler extends PixelSampler {
  final int tileSize;
  final bool randomize;

  TilePixelSampler({this.tileSize: 32, this.randomize: true});

  void setup(int x, int y, int width, int height) {
    super.setup(x, y, width, height);
    _numSamples = width * height;
    _samples = new Int32List(_numSamples * 2);
    int numXTiles = width ~/ tileSize + ((width % tileSize == 0) ? 0 : 1);
    int numYTiles = height ~/ tileSize + ((height % tileSize == 0) ? 0 : 1);

    Int32List tiles = new Int32List(numXTiles * numYTiles * 2);
    for (int yi = 0, ti = 0; yi < numYTiles; ++yi) {
      for (int xi = 0; xi < numXTiles; ++xi) {
        tiles[ti++] = xi;
        tiles[ti++] = yi;
      }
    }

    final int numTiles = tiles.length ~/ 2;

    if (randomize) {
      RNG rng = new RNG();

      // Shuffle the tiles
      for (int ti = 01; ti < numTiles; ++ti) {
        int lx = ti * 2;
        int ly = lx + 1;
        int rx = (rng.randomUint() % numTiles) * 2;
        int ry = rx + 1;

        int t = tiles[lx];
        tiles[lx] = tiles[rx];
        tiles[rx] = t;

        t = tiles[ly];
        tiles[ly] = tiles[ry];
        tiles[ry] = t;
      }
    }

    int si = 0;
    for (int i = 0, ti = 0; i < numTiles; ++i) {
      int tx = tiles[ti++];
      int ty = tiles[ti++];

      int sx = left + (tx * tileSize);
      int sy = top + (ty * tileSize);

      for (int yi = 0; yi < tileSize; ++yi) {
        int y = sy + yi;
        if (y > bottom) {
          break;
        }
        for (int xi = 0; xi < tileSize; ++xi) {
          int x = sx + xi;
          if (x > right) {
            break;
          }

          _samples[si++] = x;
          _samples[si++] = y;
        }
      }
    }
  }

  int numPixels() => _numSamples;

  void getPixel(int index, List<int> pixel) {
    index *= 2;
    if (index >= _samples.length - 1) {
      return;
    }
    pixel[0] = _samples[index];
    pixel[1] = _samples[index + 1];
  }

  static TilePixelSampler Create(ParamSet params) {
    int tileSize = params.findOneInt('tilesize', 32);
    bool randomize = params.findOneBool('randomize', true);

    return new TilePixelSampler(tileSize: tileSize, randomize: randomize);
  }

  int _numSamples;
  Int32List _samples;
}
