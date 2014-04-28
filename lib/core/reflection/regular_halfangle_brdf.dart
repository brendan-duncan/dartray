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

class RegularHalfangleBRDF extends BxDF {
  RegularHalfangleBRDF(this.brdf, this.nThetaH, this.nThetaD, this.nPhiD)
      : super(BSDF_REFLECTION | BSDF_GLOSSY);

  Spectrum f(Vector WO, Vector WI) {
    // Compute [wh] and transform wi to halfangle coordinate system
    Vector wo = new Vector.from(WO);
    Vector wi = new Vector.from(WI);
    Vector wh = wo + wi;
    if (wh.z < 0.0) {
      wo = -wo;
      wi = -wi;
      wh = -wh;
    }

    if (wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0) {
      return new Spectrum(0.0);
    }

    wh = Vector.Normalize(wh);

    double whTheta = Vector.SphericalTheta(wh);
    double whCosPhi = Vector.CosPhi(wh);
    double whSinPhi = Vector.SinPhi(wh);
    double whCosTheta = Vector.CosTheta(wh);
    double whSinTheta = Vector.SinTheta(wh);
    Vector whx = new Vector(whCosPhi * whCosTheta,
                            whSinPhi * whCosTheta,
                            -whSinTheta);
    Vector why = new Vector(-whSinPhi, whCosPhi, 0.0);
    Vector wd = new Vector(Vector.Dot(wi, whx),
                           Vector.Dot(wi, why),
                           Vector.Dot(wi, wh));

    // Compute _index_ into measured BRDF tables
    double wdTheta = Vector.SphericalTheta(wd);
    double wdPhi = Vector.SphericalPhi(wd);
    if (wdPhi > Math.PI) {
      wdPhi -= Math.PI;
    }

    // Compute indices _whThetaIndex_, _wdThetaIndex_, _wdPhiIndex_
    int REMAP(V, MAX, COUNT) => ((V / MAX).toInt() * COUNT).clamp(0, COUNT - 1);

    int whThetaIndex = REMAP(Math.sqrt(Math.max(0.0, whTheta / (Math.PI / 2.0))),
                             1.0, nThetaH);
    int wdThetaIndex = REMAP(wdTheta, Math.PI / 2.0, nThetaD);
    int wdPhiIndex = REMAP(wdPhi, Math.PI, nPhiD);

    int index = wdPhiIndex + nPhiD * (wdThetaIndex + whThetaIndex * nThetaD);

    return new Spectrum.rgb(brdf[3 * index],
                            brdf[3 * index + 1],
                            brdf[3 * index + 2]);
  }

  List<double> brdf;
  int nThetaH, nThetaD, nPhiD;
}
