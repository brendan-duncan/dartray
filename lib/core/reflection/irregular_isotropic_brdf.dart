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

class IrregIsotropicBRDFSample {
  IrregIsotropicBRDFSample(Point pp, Spectrum vv)
      : p = new Point.from(pp),
        v = new Spectrum.from(vv);

  Point p;
  Spectrum v;
}

class IrregularIsotropicBRDF extends BxDF {
  IrregularIsotropicBRDF(this.isoBRDFData)
      : super(BSDF_REFLECTION | BSDF_GLOSSY);

  Spectrum f(Vector wo, Vector wi) {
    Point m = BRDFRemap(wo, wi);
    double lastMaxDist2 = 0.001;
    while (true) {
      Spectrum v = new Spectrum(0.0);
      double sumWeights = 0.0;
      int nFound = 0;

      void proc(Point p, IrregIsotropicBRDFSample sample, double d2,
                double maxDist2) {
        double weight = Math.exp(-100.0 * d2);
        v += sample.v * weight;
        sumWeights += weight;
        ++nFound;
      }

      // Try to find enough BRDF samples around _m_ within search radius
      List<double> maxDist2 = [lastMaxDist2];
      isoBRDFData.lookup(m, proc, maxDist2);
      if (nFound > 2 || lastMaxDist2 > 1.5) {
        return v.clamp() / sumWeights;
      }

      lastMaxDist2 *= 2.0;
    }
  }

  KdTree isoBRDFData;
}
