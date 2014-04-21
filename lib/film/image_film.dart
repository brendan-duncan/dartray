/****************************************************************************
 *  Copyright (C) 2014 by Brendan Duncan.                                   *
 *                                                                          *
 *  This file is part of DartRay.                                           *
 *                                                                          *
 *  Licensed under the Apache License, Version 2.0 (the 'License');         *
 *  you may not use this file except in compliance with the License.        *
 *  You may obtain a copy of the License at                                 *
 *                                                                          *
 *  http://www.apache.org/licenses/LICENSE-2.0                              *
 *                                                                          *
 *  Unless required by applicable law or agreed to in writing, software     *
 *  distributed under the License is distributed on an 'AS IS' BASIS,       *
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 *  See the License for the specific language governing permissions and     *
 *  limitations under the License.                                          *
 *                                                                          *
 *   This project is based on PBRT v2 ; see http://www.pbrt.org             *
 *   pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.*
 ****************************************************************************/
part of film;

class ImageFilm extends Film {
  static ImageFilm Create(ParamSet params, Filter filter,
                          [Image image, PreviewCallback previewCallback]) {
    int xres = params.findOneInt('xresolution', 640);
    int yres = params.findOneInt('yresolution', 480);
    String filename = params.findOneString('filename', '');

    List<double> crop = params.findFloat('cropWindow');
    if (crop == null) {
      crop = [ 0.0, 1.0, 0.0, 1.0 ];
    }

    if (image != null) {
      xres = image.width;
      yres = image.height;
    }

    xres = (xres * RenderOverrides.GetResolutionScale()).toInt();
    yres = (yres * RenderOverrides.GetResolutionScale()).toInt();

    return new ImageFilm(xres, yres, filter, crop, filename,
                         image, previewCallback);
  }

  List<double> cropWindow;
  Image image;
  String filename;
  PreviewCallback previewCallback;
  OutputImage output;
  int samplesProcessed = 0;
  int previewCount;
  Filter filter;
  int xPixelStart;
  int yPixelStart;
  int xPixelCount;
  int yPixelCount;

  ImageFilm(int xres, int yres, this.filter, this.cropWindow, this.filename,
            [this.image, this.previewCallback]) :
    super(xres, yres),
    _gamma = new Uint8List(256) {
    double gamma = 1.0 / 2.2;
    for (int i = 0; i < 256; ++i) {
      _gamma[i] = (pow(i / 255.0, gamma) * 255.0).floor().clamp(0, 255);
    }

    // Compute film image extent
    xPixelStart = (xResolution * cropWindow[0]).ceil();
    xPixelCount = max(1, (xResolution * cropWindow[1]).ceil() - xPixelStart);
    yPixelStart = (yResolution * cropWindow[2]).ceil();
    yPixelCount = max(1, (yResolution * cropWindow[3]).ceil() - yPixelStart);

    previewCount = xPixelCount * 12;

    // Allocate film image storage
    _Lxyz = new Float32List(xPixelCount * yPixelCount * 3);
    _splatXYZ = new Float32List(xPixelCount * yPixelCount * 3);
    _weightSum = new Float32List(xPixelCount * yPixelCount);

    // Precompute filter weight table
    _filterTable = new Float32List(FILTER_TABLE_SIZE * FILTER_TABLE_SIZE);
    int fi = 0;
    for (int y = 0; y < FILTER_TABLE_SIZE; ++y) {
      double fy = (y + 0.5) * filter.yWidth / FILTER_TABLE_SIZE;
      for (int x = 0; x < FILTER_TABLE_SIZE; ++x) {
        double fx = (x + 0.5) * filter.xWidth / FILTER_TABLE_SIZE;
        _filterTable[fi++] = filter.evaluate(fx, fy);
      }
    }

    if (image == null) {
      image = new Image(xPixelCount, yPixelCount);
      image.fill(0xff888888);
    }

    output = new OutputImage(xPixelStart, yPixelStart,
                             xPixelCount, yPixelCount);
  }

  void addSample(CameraSample sample, Spectrum L) {
    // Compute sample's raster extent
    double dimageX = sample.imageX - 0.5;
    double dimageY = sample.imageY - 0.5;
    int x0 = (dimageX - filter.xWidth).ceil();
    int x1 = (dimageX + filter.xWidth).floor();
    int y0 = (dimageY - filter.yWidth).ceil();
    int y1 = (dimageY + filter.yWidth).floor();

    x0 = max(x0, xPixelStart);
    x1 = min(x1, xPixelStart + xPixelCount - 1);
    y0 = max(y0, yPixelStart);
    y1 = min(y1, yPixelStart + yPixelCount - 1);

    if ((x1 - x0) < 0 || (y1 - y0) < 0) {
      return;
    }

    // Loop over filter support and add sample to pixel arrays
    XYZColor xyz = L.toXYZ();

    // Precompute $x$ and $y$ filter table offsets
    Int32List ifx = new Int32List(x1 - x0 + 1);
    for (int x = x0; x <= x1; ++x) {
      double fx = ((x - dimageX) * filter.invXWidth * FILTER_TABLE_SIZE).abs();
      ifx[x - x0] = min(fx.floor(), FILTER_TABLE_SIZE - 1);
    }

    Int32List ify = new Int32List(y1 - y0 + 1);
    for (int y = y0; y <= y1; ++y) {
      double fy = ((y - dimageY) * filter.invYWidth * FILTER_TABLE_SIZE).abs();
      ify[y - y0] = min(fy.floor(), FILTER_TABLE_SIZE - 1);
    }

    Uint8List pixels = image.getBytes();
    List<double> rgb = [0.0, 0.0, 0.0];
    List<double> splatRGB = [0.0, 0.0, 0.0];
    for (int y = y0; y <= y1; ++y) {
      int oi = (y0 * xPixelCount + x0) * 4;
      for (int x = x0; x <= x1; ++x, oi += 4) {
        // Evaluate filter value at $(x,y)$ pixel
        int offset = ify[y - y0] * FILTER_TABLE_SIZE + ifx[x - x0];
        double filterWt = _filterTable[offset];

        // Update pixel values with filtered sample contribution
        int pi = (y - yPixelStart) * xPixelCount + (x - xPixelStart);
        int pi3 = pi * 3;

        _Lxyz[pi3] += filterWt * xyz[0];
        _Lxyz[pi3 + 1] += filterWt * xyz[1];
        _Lxyz[pi3 + 2] += filterWt * xyz[2];
        _weightSum[pi] += filterWt;

        // Update the image buffer so it can be visualized.

        // Convert pixel XYZ color to RGB
        Spectrum.XYZToRGB(_Lxyz[pi3], _Lxyz[pi3 + 1], _Lxyz[pi3 + 2], rgb);

        // Normalize pixel with weight sum
        double weightSum = _weightSum[pi];
        if (weightSum != 0.0) {
          double invWt = 1.0 / weightSum;
          rgb[0] = max(0.0, rgb[0] * invWt);
          rgb[1] = max(0.0, rgb[1] * invWt);
          rgb[2] = max(0.0, rgb[2] * invWt);
        }

        // Add splat value at pixel
        Spectrum.XYZToRGB(_splatXYZ[pi3], _splatXYZ[pi3 + 1], _splatXYZ[pi3 + 2],
                 splatRGB);

        rgb[0] += splatRGB[0];
        rgb[1] += splatRGB[1];
        rgb[2] += splatRGB[2];

        pixels[oi] = _gamma[(rgb[0] * 255.0).floor().clamp(0, 255)];
        pixels[oi + 1] = _gamma[(rgb[1] * 255.0).floor().clamp(0, 255)];
        pixels[oi + 2] = _gamma[(rgb[2] * 255.0).floor().clamp(0, 255)];
        pixels[oi + 3] = 255;
      }
    }

    samplesProcessed++;
    if (previewCallback != null && samplesProcessed % previewCount == 0) {
      previewCallback(image);
    }
  }

  void splat(CameraSample sample, Spectrum L) {
    if (L.hasNaNs()) {
      LogWarning('ImageFilm ignoring splatted spectrum with NaN values');
      return;
    }

    XYZColor xyz = L.toXYZ();

    int x = (sample.imageX).floor();
    int y = (sample.imageY).floor();

    if (x < xPixelStart || x - xPixelStart >= xPixelCount ||
        y < yPixelStart || y - yPixelStart >= yPixelCount) {
      return;
    }

    int pi = (y - yPixelStart) * xPixelCount + (x - xPixelStart);
    int pi3 = pi * 3;

    _splatXYZ[pi3] += xyz[0];
    _splatXYZ[pi3 + 1] += xyz[1];
    _splatXYZ[pi3 + 2] += xyz[2];


    // Update the image buffer so it can be visualized.
    List<double> rgb = [0.0, 0.0, 0.0];
    List<double> splatRGB = [0.0, 0.0, 0.0];
    // Convert pixel XYZ color to RGB
    Spectrum.XYZToRGB(_Lxyz[pi3], _Lxyz[pi3 + 1], _Lxyz[pi3 + 2], rgb);

    // Normalize pixel with weight sum
    double weightSum = _weightSum[pi];
    if (weightSum != 0.0) {
      double invWt = 1.0 / weightSum;
      rgb[0] = max(0.0, rgb[0] * invWt);
      rgb[1] = max(0.0, rgb[1] * invWt);
      rgb[2] = max(0.0, rgb[2] * invWt);
    }

    // Add splat value at pixel
    Spectrum.XYZToRGB(_splatXYZ[pi3], _splatXYZ[pi3 + 1], _splatXYZ[pi3 + 2],
                      splatRGB);
    rgb[0] += splatRGB[0];
    rgb[1] += splatRGB[1];
    rgb[2] += splatRGB[2];

    Uint8List pixels = image.getBytes();
    int oi = pi * 4;
    pixels[oi] = _gamma[(rgb[0] * 255.0).floor().clamp(0, 255)];
    pixels[oi + 1] = _gamma[(rgb[1] * 255.0).floor().clamp(0, 255)];
    pixels[oi + 2] = _gamma[(rgb[2] * 255.0).floor().clamp(0, 255)];
    pixels[oi + 3] = 255;

    samplesProcessed++;
    if (previewCallback != null && samplesProcessed % previewCount == 0) {
      previewCallback(image);
    }
  }

  void getSampleExtent(List<int> extent) {
    extent[0] = (xPixelStart + 0.5 - filter.xWidth).floor();
    extent[1] = (xPixelStart + 0.5 + xPixelCount + filter.xWidth).ceil();
    extent[2] = (yPixelStart + 0.5 - filter.yWidth).floor();
    extent[3] = (yPixelStart + 0.5 + yPixelCount + filter.yWidth).ceil();
  }

  void getPixelExtent(List<int> extent) {
    extent[0] = xPixelStart;
    extent[1] = xPixelStart + xPixelCount;
    extent[2] = yPixelStart;
    extent[3] = yPixelStart + yPixelCount;
  }

  void updateDisplay(int x0, int y0, int x1, int y1,
                     [double splatScale = 1.0]) {
    if (previewCallback != null) {
      previewCallback(image);
    }
  }

  OutputImage writeImage([double splatScale = 1.0]) {
    LogInfo('WRITE IMAGE $splatScale');
    // Convert image to RGB and compute final pixel values
    List<double> c = [0.0, 0.0, 0.0];
    List<double> splatRGB = [0.0, 0.0, 0.0];
    int pi = 0;
    int pi3 = 0;
    for (int y = 0; y < yPixelCount; ++y) {
      for (int x = 0; x < xPixelCount; ++x, ++pi, pi3 += 3) {
        // Convert pixel XYZ color to RGB
        Spectrum.XYZToRGB(_Lxyz[pi3], _Lxyz[pi3 + 1], _Lxyz[pi3 + 2], c);

        // Normalize pixel with weight sum
        double weightSum = _weightSum[pi];
        if (weightSum != 0.0) {
          double invWt = 1.0 / weightSum;
          output.rgb[pi3] = max(0.0, c[0] * invWt);
          output.rgb[pi3 + 1] = max(0.0, c[1] * invWt);
          output.rgb[pi3 + 2] = max(0.0, c[2] * invWt);
        }

        // Add splat value at pixel
        Spectrum.XYZToRGB(_splatXYZ[pi3], _splatXYZ[pi3 + 1],
                          _splatXYZ[pi3 + 2], splatRGB);

        output.rgb[pi3] += splatScale * splatRGB[0];
        output.rgb[pi3 + 1] += splatScale * splatRGB[1];
        output.rgb[pi3 + 2] += splatScale * splatRGB[2];
      }
    }

    return output;
  }

  Float32List _Lxyz;
  Float32List _splatXYZ;
  Float32List _weightSum;
  Float32List _filterTable;
  final Uint8List _gamma;

  static const int FILTER_TABLE_SIZE = 16;
}

