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
part of materials;

/**
  File format descriptions:

  -- Irregularly Sampled Isotropic BRDF --

  This is the format of the BRDFs in the scenes/brdfs/ folder of the pbrt
  distribution.  This is a simple text format of numbers in a particular
  format; the hash character # is used to denote a comment that continues
  to the end of the current line.

  The first number in the file gives the number of wavelengths at which the
  reflection data was measured, numWls.  This is followed by numWls values
  that give the frequency in nm of these wavelengths.  Each BRDF
  measurement is represented by 4+numWls values.  The first two give the
  (theta,phi) angles of the incident illumination direction, the next two
  give (theta,phi) for the measured reflection direction, and the following
  numWls give the spectral coefficients for the measurement at the
  wavelength specified at the start of the file.


  -- Regular Half-Angle BRDF --
  This is the file format used in the MERL BRDF database; see http://merl.com/brdf.

  This file format is a binary format, with numbers encoded in low-endian
  form.  It represents a regular 3D tabularization of BRDF samples in RGB
  color where the dimensions indexed over are (delta phi, delta theta,
  sqrt(theta_h)).  Here, theta_h is the angle between the halfangle vector
  and the normal, and delta theta and delta phi are the offset in theta and
  phi of one of the two directions.  (Note that the offset would be the
  same for the other direction, since it's from the half-angle vector.)

  The starts with three 32-bit integers, giving the resolution of the
  overall table.  It then containes a number of samples equal to the
  product of those three integers, times 3 (for RGB).  Samples are laid out
  with delta phi the minor index, then delta theta, then sqrt(theta_h) as
  the major index.

  In the file each sample should be scaled by RGB(1500,1500,1500/1.6) of
  the original measurement.  (In order words, the sample values are scaled
  by the inverse of that as they are read in.
*/
class MeasuredMaterial extends Material {
  static MeasuredMaterial Create(Transform xform, TextureParams mp) {
    Texture bumpMap = mp.getFloatTextureOrNull('bumpmap');
    return new MeasuredMaterial(mp.findFilename('filename'), bumpMap);
  }

  MeasuredMaterial(String filename, this.bumpMap) {
    String suffix = filename.substring(filename.lastIndexOf('.')).toLowerCase();

    if (suffix.isEmpty) {
      LogError('No suffix in measured BRDF filename "$filename". '
               'Can\'t determine file type (.brdf / .merl)');
      return;
    }

    if (suffix == '.brdf') {
      // Load $(\theta, \phi)$ measured BRDF data
      if (loadedThetaPhi.containsKey(filename)) {
        thetaPhiData = loadedThetaPhi[filename];
        return;
      }

      Completer c = new Completer();
      ResourceManager.RequestFile(filename, c.future).then((bytes) {
        List<double> values = ReadFloatFile(bytes, filename);

        int pos = 0;
        int numWls = values[pos++].toInt();
        if ((values.length - 1 - numWls) % (4 + numWls) != 0) {
          LogError('Excess or insufficient data in theta, phi '
                   'BRDF file \"$filename\"');
          c.complete();
          return;
        }

        Float32List wls = new Float32List(numWls);
        for (int i = 0; i < numWls; ++i) {
          wls[i] = values[pos++];
        }

        BBox bbox = new BBox();
        List<IrregIsotropicBRDFSample> samples = [];
        while (pos < values.length) {
          double thetai = values[pos++];
          double phii = values[pos++];
          double thetao = values[pos++];
          double phio = values[pos++];
          Vector wo = Vector.SphericalDirection(sin(thetao), cos(thetao), phio);
          Vector wi = Vector.SphericalDirection(sin(thetai), cos(thetai), phii);
          Spectrum s = new Spectrum.fromSampled(wls, values, pos);

          pos += numWls;

          Point p = BRDFRemap(wo, wi);
          samples.add(new IrregIsotropicBRDFSample(p, s));
          bbox = BBox.UnionPoint(bbox, p);
        }

        loadedThetaPhi[filename] = thetaPhiData = new KdTree(samples);
        c.complete();
      });
    } else {
      // Load RegularHalfangle BRDF Data
      nThetaH = 90;
      nThetaD = 90;
      nPhiD = 180;

      if (loadedRegularHalfangle.containsKey(filename)) {
        regularHalfangleData = loadedRegularHalfangle[filename];
        return;
      }

      Completer c = new Completer();
      ResourceManager.RequestFile(filename, c.future).then((bytes) {
        InputBuffer fp = new InputBuffer(bytes, bigEndian: false);
        int dims0 = fp.readInt32();
        int dims1 = fp.readInt32();
        int dims2 = fp.readInt32();

        int n = dims0 * dims1 * dims2;
        if (n != nThetaH * nThetaD * nPhiD) {
          LogError('Dimensions don\'t match');
          c.complete();
          return;
        }

        regularHalfangleData = new Float32List(3 * n);
        final int chunkSize = 2 * nPhiD;
        Float32List tmp = new Float32List(chunkSize);
        int nChunks = n ~/ chunkSize;
        assert((n % chunkSize) == 0);

        const List<double> scales = const [1.0 / 1500.0,
                                           1.15 / 1500.0,
                                           1.66 / 1500.0];

        // There isn't a readFloat64 method in InputBuffer, so read the
        // bits as a Uint64, and typecast it to a Float64 (double) using
        // a TypedData.view.
        Uint64List i64 = new Uint64List(1);
        Float64List f64 = new Float64List.view(i64.buffer);

        for (int c = 0; c < 3; ++c) {
          int offset = 0;
          for (int i = 0; i < nChunks; ++i) {
            for (int j = 0; j < chunkSize; ++j) {
              i64[0] = fp.readUint64();
              tmp[j] = f64[0];
            }

            for (int j = 0; j < chunkSize; ++j) {
              double t = tmp[j] * scales[c];
              regularHalfangleData[3 * offset + c] = max(0.0, t);
              offset++;
            }
          }
        }

        loadedRegularHalfangle[filename] = regularHalfangleData;
      });
    }
  }

  BSDF getBSDF(DifferentialGeometry dgGeom, DifferentialGeometry dgShading) {
    DifferentialGeometry dgs = new DifferentialGeometry();
    if (bumpMap != null) {
      Material.Bump(bumpMap, dgGeom, dgShading, dgs);
    } else {
      dgs = dgShading;
    }

    BSDF bsdf = new BSDF(dgs, dgGeom.nn);

    if (regularHalfangleData != null) {
      bsdf.add(new RegularHalfangleBRDF(regularHalfangleData, nThetaH,
                                        nThetaD, nPhiD));
    } else if (thetaPhiData != null) {
      bsdf.add(new IrregularIsotropicBRDF(thetaPhiData));
    }

    return bsdf;
  }

  KdTree thetaPhiData;
  List<double> regularHalfangleData;
  int nThetaH, nThetaD, nPhiD;
  Texture bumpMap;

  static Map<String, Float32List> loadedRegularHalfangle = {};
  static Map<String, KdTree> loadedThetaPhi = {};
}
