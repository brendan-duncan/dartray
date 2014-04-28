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

const int BSDF_REFLECTION = 1 << 0;
const int BSDF_TRANSMISSION = 1<<1;
const int BSDF_DIFFUSE = 1 << 2;
const int BSDF_GLOSSY = 1 << 3;
const int BSDF_SPECULAR = 1 << 4;
const int BSDF_ALL_TYPES = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR;
const int BSDF_ALL_REFLECTION = BSDF_REFLECTION | BSDF_ALL_TYPES;
const int BSDF_ALL_TRANSMISSION = BSDF_TRANSMISSION | BSDF_ALL_TYPES;
const int BSDF_ALL = BSDF_ALL_REFLECTION | BSDF_ALL_TRANSMISSION;

/**
 * A Bidirectional Scattering Distribution Function ([BSDF]) describes the
 * way light is scattered by a surface.
 *
 * It is defined by a set of [BxDF] functions, including
 * Bidirectional Reflectance Distribution Functions [BRDF], and
 * Bidirectional Transmission Distribution Functions [BTDF].
 */
class BSDF {
  DifferentialGeometry dgShading;
  double eta;

  BSDF(this.dgShading, Normal ngeom, [this.eta = 1.0]) {
    ng = ngeom;
    nn = dgShading.nn;
    sn = Vector.Normalize(dgShading.dpdu);
    tn = Vector.Cross(nn, sn);
    nBxDFs = 0;
  }

  Spectrum sample_f(Vector woW, Vector wiW, BSDFSample bsdfSample,
                    List<double> pdf, [int flags = BSDF_ALL,
                    List<int> sampledType]) {
    Stats.STARTED_BSDF_SAMPLE();
    // Choose which _BxDF_ to sample
    int matchingComps = numComponents(flags);
    if (matchingComps == 0) {
      pdf[0] = 0.0;
      if (sampledType != null) {
        sampledType[0] = 0;
      }
      Stats.FINISHED_BSDF_SAMPLE();
      return new RGBColor(0.0);
    }

    int which = Math.min((bsdfSample.uComponent * matchingComps).floor(),
                         matchingComps - 1);

    BxDF bxdf;
    int count = which;
    for (int i = 0; i < nBxDFs; ++i) {
      if (bxdfs[i].matchesFlags(flags) && count-- == 0) {
        bxdf = bxdfs[i];
        break;
      }
    }
    assert(bxdf != null);

    // Sample chosen _BxDF_
    Vector wo = worldToLocal(woW);
    Vector wi = new Vector();
    pdf[0] = 0.0;

    Spectrum f = bxdf.sample_f(wo, wi, bsdfSample.uDir[0],
                               bsdfSample.uDir[1], pdf);
    if (pdf[0] == 0.0) {
      if (sampledType != null) {
        sampledType[0] = 0;
      }
      Stats.FINISHED_BSDF_SAMPLE();
      return new Spectrum(0.0);
    }

    if (sampledType != null) {
      sampledType[0] = bxdf.type;
    }

    wiW.copy(localToWorld(wi));

    // Compute overall PDF with all matching _BxDF_s
    if (!(bxdf.type & BSDF_SPECULAR != 0) && matchingComps > 1) {
      for (int i = 0; i < nBxDFs; ++i) {
        if (bxdfs[i] != bxdf && bxdfs[i].matchesFlags(flags)) {
          pdf[0] += bxdfs[i].pdf(wo, wi);
        }
      }
    }

    if (matchingComps > 1) {
      pdf[0] /= matchingComps;
    }

    // Compute value of BSDF for sampled direction
    if ((bxdf.type & BSDF_SPECULAR) == 0) {
      f = new Spectrum(0.0);
      if (Vector.Dot(wiW, ng) * Vector.Dot(woW, ng) > 0) {
        flags = flags & ~BSDF_TRANSMISSION;
      } else { // ignore BRDFs
        flags = flags & ~BSDF_REFLECTION;
      }

      for (int i = 0; i < nBxDFs; ++i) {
        if (bxdfs[i].matchesFlags(flags)) {
          f += bxdfs[i].f(wo, wi);
        }
      }
    }

    Stats.FINISHED_BSDF_SAMPLE();
    return f;
  }

  double pdf(Vector woW, Vector wiW, [int flags = BSDF_ALL]) {
    if (nBxDFs == 0.0) {
      return 0.0;
    }
    Stats.STARTED_BSDF_PDF();

    Vector wo = worldToLocal(woW);
    Vector wi = worldToLocal(wiW);
    double pdf = 0.0;
    int matchingComps = 0;
    for (int i = 0; i < nBxDFs; ++i) {
      if (bxdfs[i].matchesFlags(flags)) {
        ++matchingComps;
        pdf += bxdfs[i].pdf(wo, wi);
      }
    }

    double v = matchingComps > 0 ? pdf / matchingComps : 0.0;
    Stats.FINISHED_BSDF_PDF();

    return v;
  }

  void add(BxDF bxdf) {
    bxdfs[nBxDFs++] = bxdf;
  }

  int numComponents([int flags]) {
    if (flags == null) {
      return nBxDFs;
    }

    int num = 0;
    for (int i = 0; i < nBxDFs; ++i) {
      if (bxdfs[i].matchesFlags(flags)) {
        ++num;
      }
    }

    return num;
  }

  Vector worldToLocal(Vector v) {
    return new Vector(Vector.Dot(v, sn), Vector.Dot(v, tn), Vector.Dot(v, nn));
  }

  Vector localToWorld(Vector v) {
    return new Vector(sn.x * v.x + tn.x * v.y + nn.x * v.z,
                      sn.y * v.x + tn.y * v.y + nn.y * v.z,
                      sn.z * v.x + tn.z * v.y + nn.z * v.z);
  }

  Spectrum f(Vector woW, Vector wiW, [int flags = BSDF_ALL]) {
    Stats.STARTED_BSDF_EVAL();
    Vector wi = worldToLocal(wiW);
    Vector wo = worldToLocal(woW);

    if (Vector.Dot(wiW, ng) * Vector.Dot(woW, ng) > 0) {
      // ignore BTDFs
      flags = flags & ~BSDF_TRANSMISSION;
    } else {
      // ignore BRDFs
      flags = flags & ~BSDF_REFLECTION;
    }

    Spectrum f = new Spectrum(0.0);

    for (int i = 0; i < nBxDFs; ++i) {
      if (bxdfs[i].matchesFlags(flags)) {
        f += bxdfs[i].f(wo, wi);
      }
    }

    Stats.FINISHED_BSDF_EVAL();

    return f;
  }

  Spectrum rho(RNG rng, [int flags = BSDF_ALL, int sqrtSamples = 6]) {
    int nSamples = sqrtSamples * sqrtSamples;

    List<double> s1 = new List<double>(2 * nSamples);
    StratifiedSample2D(s1, sqrtSamples, sqrtSamples, rng);

    List<double> s2 = new List<double>(2 * nSamples);
    StratifiedSample2D(s2, sqrtSamples, sqrtSamples, rng);

    Spectrum ret = new Spectrum(0.0);
    for (int i = 0; i < nBxDFs; ++i) {
      if (bxdfs[i].matchesFlags(flags)) {
        ret += bxdfs[i].rho2(nSamples, s1, s2);
      }
    }

    return ret;
  }

  Spectrum rho2(Vector wo, RNG rng, [int flags = BSDF_ALL,
                int sqrtSamples = 6]) {
    int nSamples = sqrtSamples * sqrtSamples;
    List<double> s1 = new List<double>(2 * nSamples);

    StratifiedSample2D(s1, sqrtSamples, sqrtSamples, rng);

    Spectrum ret = new Spectrum(0.0);
    for (int i = 0; i < nBxDFs; ++i) {
      if (bxdfs[i].matchesFlags(flags)) {
        ret += bxdfs[i].rho(wo, nSamples, s1);
      }
    }

    return ret;
  }

  Normal nn, ng;
  Vector sn, tn;
  int nBxDFs;

  static const int MAX_BxDFS = 8;
  List<BxDF> bxdfs = new List<BxDF>(MAX_BxDFS);
}
