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
 * Abstract base class for defining both BRDF reflectance functions and
 * BTDF transmission functions. The type flags are used to identify whether
 * the BxDF is a BRDF or a BTDF.
 */
abstract class BxDF {
  BxDF(this.type);

  bool matchesFlags(int flags) {
    return (type & flags) == type;
  }

  Spectrum f(Vector wo, Vector wi);

  Spectrum sample_f(Vector wo, Vector out_wi,
                    double u1, double u2, List<double> out_pdf) {
    // Cosine-sample the hemisphere, flipping the direction if necessary
    out_wi.copy(CosineSampleHemisphere(u1, u2));
    if (wo.z < 0.0) {
      out_wi.z *= -1.0;
    }

    out_pdf[0] = pdf(wo, out_wi);

    return f(wo, out_wi);
  }

  Spectrum rho(Vector w, int nSamples, List<double> samples) {
    Spectrum r = new Spectrum(0.0);
    for (int i = 0; i < nSamples; ++i) {
      // Estimate one term of $\rho_\roman{hd}$
      Vector wi = new Vector();
      List<double> pdf = [0.0];
      Spectrum f = sample_f(w, wi, samples[2 * i], samples[2 * i + 1], pdf);
      if (pdf[0] > 0.0) {
        r += f * (Vector.AbsCosTheta(wi) / pdf[0]);
      }
    }
    return r / nSamples;
  }

  Spectrum rho2(int nSamples, List<double> samples1, List<double> samples2) {
    Spectrum r = new Spectrum(0.0);
    for (int i = 0; i < nSamples; ++i) {
      Vector wi = new Vector();
      Vector wo = UniformSampleHemisphere(samples1[2 * i], samples1[2 * i + 1]);

      double pdf_o = INV_TWOPI;
      List<double> pdf_i = [0.0];

      Spectrum f = sample_f(wo, wi, samples2[2 * i], samples2[2 * i + 1],
                            pdf_i);
      if (pdf_i[0] > 0.0) {
        r += f * (Vector.AbsCosTheta(wi) * Vector.AbsCosTheta(wo) /
            (pdf_o * pdf_i[0]));
      }
    }

    return r / (Math.PI * nSamples);
  }

  double pdf(Vector wo, Vector wi) {
    return Vector.SameHemisphere(wo, wi) ?
           Vector.AbsCosTheta(wi) * INV_PI :
           0.0;
  }

  final int type;
}
