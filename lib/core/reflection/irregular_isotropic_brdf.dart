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

class IrregIsotropicBRDFSample {
  IrregIsotropicBRDFSample(Point pp, RGBColor vv) :
    p = new Point.from(pp),
    v = new RGBColor.from(vv);

  Point p;
  RGBColor v;
}

class IrregularIsotropicBRDF extends BxDF {
  IrregularIsotropicBRDF(this.isoBRDFData) :
    super(BSDF_REFLECTION | BSDF_GLOSSY);

  RGBColor f(Vector wo, Vector wi) {
    Point m = _BRDFRemap(wo, wi);
    double lastMaxDist2 = 0.001;
    while (true) {
      RGBColor v = new RGBColor(0.0);
      double sumWeights = 0.0;
      int nFound = 0;

      void proc(Point p, IrregIsotropicBRDFSample sample, double d2,
                double maxDist2) {
        double weight = Math.exp(-100.0 * d2);
        v += sample.v.scaled(weight);
        sumWeights += weight;
        ++nFound;
      }

      // Try to find enough BRDF samples around _m_ within search radius
      double maxDist2 = lastMaxDist2;
      isoBRDFData.lookup(m, proc, maxDist2);
      if (nFound > 2 || lastMaxDist2 > 1.5) {
        return v.clamp() / sumWeights;
      }

      lastMaxDist2 *= 2.0;
    }
  }

  Point _BRDFRemap(Vector wo, Vector wi) {
    double cosi = Vector.CosTheta(wi);
    double coso = Vector.CosTheta(wo);
    double sini = Vector.SinTheta(wi);
    double sino = Vector.SinTheta(wo);
    double phii = Vector.SphericalPhi(wi);
    double phio = Vector.SphericalPhi(wo);
    double dphi = phii - phio;

    if (dphi < 0.0) {
      dphi += 2.0 * Math.PI;
    }

    if (dphi > 2.0 * Math.PI) {
      dphi -= 2.0 * Math.PI;
    }

    if (dphi > Math.PI) {
      dphi = 2.0 * Math.PI - dphi;
    }

    return new Point(sini * sino, dphi / Math.PI, cosi * coso);
  }

  KdTree isoBRDFData;
}
