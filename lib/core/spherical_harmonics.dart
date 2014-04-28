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
 * A namespace container for functions related to spherical harmonics.
 */
class SphericalHarmonics {
  static int Terms(int lmax) => (lmax + 1) * (lmax + 1);

  static int Index(int l, int m) =>
    l * l + l + m;

  static void Evaluate(Vector w, int lmax, List<double> out,
                       [int outIndex = 0]) {
    if (lmax > 28) {
        LogSevere('SHEvaluate() runs out of numerical precision for lmax > 28.'
                  'If you need more bands, try recompiling using doubles.');
    }

    // Compute Legendre polynomial values for $\cos\theta$
    assert(w.length() > 0.995 && w.length() < 1.005);
    _legendrep(w.z, lmax, out, outIndex);

    // Compute $K_l^m$ coefficients
    List<double> Klm = new List<double>(Terms(lmax));

    for (int l = 0; l <= lmax; ++l) {
      for (int m = -l; m <= l; ++m) {
        Klm[Index(l, m)] = _K(l, m);
      }
    }

    // Compute $\sin\phi$ and $\cos\phi$ values
    List<double> sins = new List<double>(lmax + 1);
    List<double> coss = new List<double>(lmax + 1);
    double xyLen = Math.sqrt(Math.max(0.0, 1.0 - w.z * w.z));
    if (xyLen == 0.0) {
      for (int i = 0; i <= lmax; ++i) {
        sins[i] = 0.0;
      }
      for (int i = 0; i <= lmax; ++i) {
        coss[i] = 1.0;
      }
    } else {
      _sinCosIndexed(w.y / xyLen, w.x / xyLen, lmax + 1, sins, coss);
    }

    // Apply SH definitions to compute final $(l,m)$ values
    final double sqrt2 = Math.sqrt(2.0);

    for (int l = 0; l <= lmax; ++l) {
      for (int m = -l; m < 0; ++m) {
        out[Index(l, m)] = sqrt2 * Klm[Index(l, m)] *
            out[Index(l, -m)] * sins[-m];
        assert(out[Index(l,m)].isNaN);
        assert(out[Index(l,m)].isFinite);
      }

      out[Index(l, 0)] *= Klm[Index(l, 0)];
      for (int m = 1; m <= l; ++m) {
        out[Index(l, m)] *= sqrt2 * Klm[Index(l, m)] * coss[m];
        assert(!out[Index(l, m)].isNaN);
        assert(!out[Index(l, m)].isFinite);
      }
    }
  }

  static void ProjectCube(func, Point p, int res, int lmax,
                          List<Spectrum> coeffs) {
    List<double> Ylm = new List<double>(Terms(lmax));

    for (int u = 0; u < res; ++u) {
      double fu = -1.0 + 2.0 * (u + 0.5) / res;
      for (int v = 0; v < res; ++v) {
        double fv = -1.0 + 2.0 * (v + 0.5) / res;
        // Incorporate results from $+z$ face to coefficients
        Vector w = new Vector(fu, fv, 1.0);
        Evaluate(Vector.Normalize(w), lmax, Ylm);
        Spectrum f = func(u, v, p, w);
        double dA = 1.0 / Math.pow(Vector.Dot(w, w), 3.0 / 2.0);
        for (int k = 0; k < Terms(lmax); ++k) {
          coeffs[k] += f * dA * Ylm[k] * (4.0 / (res * res));
        }

        // Incorporate results from other faces to coefficients
        w = new Vector(fu, fv, -1.0);
        Evaluate(Vector.Normalize(w), lmax, Ylm);
        f = func(u, v, p, w);
        for (int k = 0; k < Terms(lmax); ++k) {
          coeffs[k] += f * Ylm[k] * dA * (4.0 / (res * res));
        }

        w = new Vector(fu, 1.0, fv);
        Evaluate(Vector.Normalize(w), lmax, Ylm);
        f = func(u, v, p, w);
        for (int k = 0; k < Terms(lmax); ++k) {
          coeffs[k] += f * Ylm[k] * dA * (4.0 / (res * res));
        }

        w = new Vector(fu, -1.0, fv);
        Evaluate(Vector.Normalize(w), lmax, Ylm);
        f = func(u, v, p, w);
        for (int k = 0; k < Terms(lmax); ++k) {
          coeffs[k] += f * Ylm[k] * dA * (4.0 / (res * res));
        }

        w = new Vector(1.0, fu, fv);
        Evaluate(Vector.Normalize(w), lmax, Ylm);
        f = func(u, v, p, w);
        for (int k = 0; k < Terms(lmax); ++k) {
          coeffs[k] += f * Ylm[k] * dA * (4.0 / (res * res));
        }

        w = new Vector(-1.0, fu, fv);
        Evaluate(Vector.Normalize(w), lmax, Ylm);
        f = func(u, v, p, w);
        for (int k = 0; k < Terms(lmax); ++k) {
          coeffs[k] += f * Ylm[k] * dA * (4.0 / (res * res));
        }
      }
    }
  }

  static void ProjectIncidentDirectRadiance(Point p, double pEpsilon,
                                             double time, Scene scene,
                                             bool computeLightVisibility,
                                             int lmax, RNG rng,
                                             List<Spectrum> c_d) {
    // Loop over light sources and sum their SH coefficients
    List<Spectrum> c = new List<Spectrum>(Terms(lmax));
    for (int i = 0, len = c.length; i < len; ++i) {
      c[i] = new Spectrum(0.0);
    }

    for (int i = 0; i < scene.lights.length; ++i) {
      Light light = scene.lights[i];
      light.shProject(p, pEpsilon, lmax, scene, computeLightVisibility, time,
                      rng, c);
       for (int j = 0; j < Terms(lmax); ++j) {
         c_d[j] += c[j];
       }
    }

    ReduceRinging(c_d, lmax);
  }

  static void ProjectIncidentIndirectRadiance(Point p, double pEpsilon,
                                              double time, Renderer renderer,
                                              Sample origSample,
                                              Scene scene, int lmax, RNG rng,
                                              int nSamples,
                                              List<Spectrum> c_i) {

    Sample sample = new Sample.from(origSample);
    List<int> scramble = [rng.randomUint(), rng.randomUint()];
    nSamples = RoundUpPow2(nSamples);
    List<double> Ylm = new List<double>(Terms(lmax));
    for (int i = 0; i < nSamples; ++i) {
      // Sample incident direction for radiance probe
      List<double> u = [0.0, 0.0];
      Sample02(i, scramble, u);
      Vector wi = UniformSampleSphere(u[0], u[1]);
      double pdf = UniformSpherePdf();

      // Compute incident radiance along direction for probe
      Spectrum Li = new Spectrum(0.0);
      RayDifferential ray = new RayDifferential(p, wi, pEpsilon, INFINITY, time);

      // Fill in values in _sample_ for radiance probe ray
      sample.time = time;
      for (int j = 0; j < sample.n1D.length; ++j) {
        for (int k = 0; k < sample.n1D[j]; ++k) {
          sample.oneD[j][k] = rng.randomFloat();
        }
      }

      for (int j = 0; j < sample.n2D.length; ++j) {
        for (int k = 0; k < 2 * sample.n2D[j]; ++k) {
          sample.twoD[j][k] = rng.randomFloat();
        }
      }

      Li = renderer.Li(scene, ray, sample, rng);

      // Update SH coefficients for probe sample point
      Evaluate(wi, lmax, Ylm);
      for (int j = 0; j < Terms(lmax); ++j) {
        c_i[j] += Li * Ylm[j] / (pdf * nSamples);
      }
    }
  }

  static void ReduceRinging(List<Spectrum> c, int lmax, [double lambda = 0.005]) {
    for (int l = 0; l <= lmax; ++l) {
      double scale = 1.0 / (1.0 + lambda * l * l * (l + 1) * (l + 1));
      for (int m = -l; m <= l; ++m) {
        c[Index(l, m)].scale(scale);
      }
    }
  }

  static void Rotate(List<Spectrum> c_in, List<Spectrum> c_out, Matrix4x4 m,
                     int lmax) {
    List<double> alpha = [0.0],
                 beta = [0.0],
                 gamma = [0.0];

    _toZYZ(m, alpha, beta, gamma);
    List<Spectrum> work = Spectrum.AllocateList(Terms(lmax));

    RotateZ(c_in, c_out, gamma[0], lmax);
    RotateXPlus(c_out, work, lmax);
    RotateZ(work, c_out, beta[0], lmax);
    RotateXMinus(c_out, work, lmax);
    RotateZ(work, c_out, alpha[0], lmax);
  }

  static void RotateZ(List<Spectrum> c_in, List<Spectrum> c_out, double alpha,
                      int lmax) {
    assert(c_in != c_out);
    c_out[0].copy(c_in[0]);
    if (lmax == 0) {
      return;
    }

    // Precompute sine and cosine terms for $z$-axis SH rotation
    List<double> ct = new List<double>(lmax + 1);
    List<double> st = new List<double>(lmax + 1);

    _sinCosIndexed(Math.sin(alpha), Math.cos(alpha), lmax + 1, st, ct);

    for (int l = 1; l <= lmax; ++l) {
      // Rotate coefficients for band _l_ about $z$
      for (int m = -l; m < 0; ++m) {
        c_out[Index(l, m)] = c_in[Index(l,  m)] * ct[-m] +
                               c_in[Index(l, -m)] * -st[-m];
      }

      c_out[Index(l, 0)].copy(c_in[Index(l, 0)]);
      for (int m = 1; m <= l; ++m) {
        c_out[Index(l, m)] = c_in[Index(l,  m)] * ct[m] +
                               c_in[Index(l, -m)] * st[m];
      }
    }
  }

  static void RotateXMinus(List<Spectrum> c_in, List<Spectrum> c_out,
                           int lmax) {
    // -x rotations are the same as +x rotations, just with a negation
    // factor thrown in for some of the terms.
    RotateXPlus(c_in, c_out, lmax);

    // l = 0 band is a no op...
    for (int l = 1; l <= lmax; ++l) {
      double s = (l & 0x1) != 0 ? -1.0 : 1.0;
      c_out[Index(l, 0)].scale(s);
      for (int m = 1; m <= l; ++m) {
        s = -s;
        c_out[Index(l, m)].scale(s);
        c_out[Index(l, -m)].scale(-s);
      }
    }
  }

  static void RotateXPlus(List<Spectrum> c_in, List<Spectrum> c_out, int lmax) {
    O(l, m) => c_in[Index(l, m)];
    int oi = 0;

    // first band is a no-op
    c_out[oi++] = c_in[0];

    if (lmax < 1) {
      return;
    }
    c_out[oi++] = O(1, 0);
    c_out[oi++] = -O(1, -1);
    c_out[oi++] = O(1, 1);

    if (lmax < 2) {
      return;
    }
    c_out[oi++] = O(2, 1);
    c_out[oi++] = -O(2, -1);
    c_out[oi++] = O(2, 0) * -0.5 + O(2, 2) * -0.8660254037844386;
    c_out[oi++] = (-O(2, -2));
    c_out[oi++] = (O(2, 0) * -0.8660254037844386 + O(2, 2) * 0.5);

    // Remainder of SH x+ rotation definition
    if (lmax < 3) {
      return;
    }
    c_out[oi++] = (O(3, 0) * -0.7905694150420949 +
                   O(3, 2) * 0.6123724356957945);
    c_out[oi++] = -O(3,-2);
    c_out[oi++] = (O(3, 0) * -0.6123724356957945 +
                   O(3,2) * -0.7905694150420949);
    c_out[oi++] = (O(3, -3) * 0.7905694150420949 +
                   O(3, -1) * 0.6123724356957945);
    c_out[oi++] = (O(3, 1) * -0.25 + O(3, 3) * -0.9682458365518543);
    c_out[oi++] = (O(3, -3) * -0.6123724356957945 +
                   O(3, -1) * 0.7905694150420949);
    c_out[oi++] = O(3, 1) * -0.9682458365518543 + O(3, 3) * 0.25;

    if (lmax < 4) {
      return;
    }
    c_out[oi++] = O(4, 1) * -0.9354143466934853 + O(4, 3) * 0.35355339059327373;
    c_out[oi++] = O(4, -3) * -0.75 + O(4,-1) * 0.6614378277661477;
    c_out[oi++] = O(4, 1) * -0.35355339059327373 +
                  O(4, 3) * -0.9354143466934853;
    c_out[oi++] = O(4, -3) * 0.6614378277661477 + O(4, -1) * 0.75;
    c_out[oi++] = O(4, 0) * 0.375 + O(4, 2) * 0.5590169943749475 +
                  O(4, 4) * 0.739509972887452;
    c_out[oi++] = O(4, -4) * 0.9354143466934853 +
                  O(4, -2) * 0.35355339059327373;
    c_out[oi++] = O(4, 0) * 0.5590169943749475 + O(4, 2) * 0.5 +
                  O(4, 4) * -0.6614378277661477;
    c_out[oi++] = O(4, -4) * -0.35355339059327373 +
                  O(4, -2) * 0.9354143466934853;
    c_out[oi++] = O(4, 0) * 0.739509972887452 + O(4, 2) * -0.6614378277661477 +
                  O(4, 4) * 0.125;

    if (lmax < 5) {
      return;
    }

    c_out[oi++] = (O(5,0) * 0.701560760020114 - O(5,2) * 0.6846531968814576 +
                O(5,4) * 0.19764235376052372);
    c_out[oi++] = (O(5,-4) * -0.5 + O(5,-2) * 0.8660254037844386);
    c_out[oi++] = (O(5,0) * 0.5229125165837972 + O(5,2) * 0.30618621784789724 -
                O(5,4) * 0.795495128834866);
    c_out[oi++] = (O(5,-4) * 0.8660254037844386 + O(5,-2) * 0.5);
    c_out[oi++] = (O(5,0) * 0.4841229182759271 + O(5,2) * 0.6614378277661477 +
                O(5,4) * 0.57282196186948);
    c_out[oi++] = (O(5,-5) * -0.701560760020114 - O(5,-3) * 0.5229125165837972 -
                O(5,-1) * 0.4841229182759271);
    c_out[oi++] = (O(5,1) * 0.125 + O(5,3) * 0.4050462936504913 + O(5,5) *
                   0.9057110466368399);
    c_out[oi++] = (O(5,-5) * 0.6846531968814576 - O(5,-3) *
                   0.30618621784789724 - O(5,-1) * 0.6614378277661477);
    c_out[oi++] = (O(5,1) * 0.4050462936504913 + O(5,3) * 0.8125 - O(5,5) *
                   0.4192627457812106);
    c_out[oi++] = (O(5,-5) * -0.19764235376052372 + O(5,-3) *
                   0.795495128834866 - O(5,-1) * 0.57282196186948);
    c_out[oi++] = (O(5,1) * 0.9057110466368399 - O(5,3) * 0.4192627457812106 +
                   O(5,5) * 0.0625);

    if (lmax < 6) {
      return;
    }
    c_out[oi++] = (O(6,1) * 0.879452954966893 - O(6,3) * 0.46351240544347894 +
                   O(6,5) * 0.10825317547305482);
    c_out[oi++] = (O(6,-5) * -0.3125 + O(6,-3) * 0.8028270361665706 -
                   O(6,-1) * 0.5077524002897476);
    c_out[oi++] = (O(6,1) * 0.4330127018922193 + O(6,3) * 0.6846531968814576 -
                O(6,5) * 0.5863019699779287);
    c_out[oi++] = (O(6,-5) * 0.8028270361665706 - O(6,-3) * 0.0625 -
                   O(6,-1) * 0.5929270612815711);
    c_out[oi++] = (O(6,1) * 0.19764235376052372 + O(6,3) * 0.5625 +
                   O(6,5) * 0.8028270361665706);
    c_out[oi++] = (O(6,-5) * -0.5077524002897476 -
                   O(6,-3) * 0.5929270612815711 - O(6,-1) * 0.625);
    c_out[oi++] = (O(6,0) * -0.3125 - O(6,2) * 0.45285552331841994 -
                   O(6,4) * 0.49607837082461076 - O(6,6) * 0.6716932893813962);
    c_out[oi++] = (O(6,-6) * -0.879452954966893 - O(6,-4) * 0.4330127018922193 -
                   O(6,-2) * 0.19764235376052372);
    c_out[oi++] = (O(6,0) * -0.45285552331841994 - O(6,2) * 0.53125 -
                   O(6,4) * 0.1711632992203644 + O(6,6) * 0.6952686081652184);
    c_out[oi++] = (O(6,-6) * 0.46351240544347894 - O(6,-4) *
                   0.6846531968814576 - O(6,-2) * 0.5625);
    c_out[oi++] = (O(6,0) * -0.49607837082461076 - O(6,2) * 0.1711632992203644 +
                   O(6,4) * 0.8125 - O(6,6) * 0.2538762001448738);
    c_out[oi++] = (O(6,-6) * -0.10825317547305482 +
                   O(6,-4) * 0.5863019699779287 - O(6,-2) * 0.8028270361665706);
    c_out[oi++] = (O(6,0) * -0.6716932893813962 + O(6,2) * 0.6952686081652184 -
                   O(6,4) * 0.2538762001448738 + O(6,6) * 0.03125);

    if (lmax < 7) {
      return;
    }
    c_out[oi++] = (O(7,0) * -0.6472598492877494 + O(7,2) * 0.6991205412874092 -
                   O(7,4) * 0.2981060004427955 + O(7,6) * 0.05846339666834283);
    c_out[oi++] = (O(7,-6) * -0.1875 + O(7,-4) * 0.6373774391990981 -
                   O(7,-2) * 0.7473912964438374);
    c_out[oi++] = (O(7,0) * -0.47495887979908324 - O(7,2) * 0.07328774624724109 +
                   O(7,4) * 0.78125 - O(7,6) * 0.3983608994994363);
    c_out[oi++] = (O(7,-6) * 0.6373774391990981 - O(7,-4) * 0.5 -
                   O(7,-2) * 0.5863019699779287);
    c_out[oi++] = (O(7,0) * -0.42961647140211 - O(7,2) * 0.41984465132951254 +
                   O(7,4) * 0.10364452469860624 + O(7,6) * 0.7927281808728639);
    c_out[oi++] = (O(7,-6) * -0.7473912964438374 - O(7,-4) * 0.5863019699779287 -
                   O(7,-2) * 0.3125);
    c_out[oi++] = (O(7,0) * -0.41339864235384227 - O(7,2) * 0.5740991584648073 -
                   O(7,4) * 0.5385527481129402 - O(7,6) * 0.4576818286211503);
    c_out[oi++] = (O(7,-7) * 0.6472598492877494 +
                   O(7,-5) * 0.47495887979908324 +
                   O(7,-3) * 0.42961647140211 + O(7,-1) * 0.41339864235384227);
    c_out[oi++] = (O(7,1) * -0.078125 - O(7,3) * 0.24356964481437335 -
                   O(7,5) * 0.4487939567607835 - O(7,7) * 0.8562442974262661);
    c_out[oi++] = (O(7,-7) * -0.6991205412874092 +
                   O(7,-5) * 0.07328774624724109 +
                   O(7,-3) * 0.41984465132951254 +
                   O(7,-1) * 0.5740991584648073);
    c_out[oi++] = (O(7,1) * -0.24356964481437335 - O(7,3) * 0.609375 -
                   O(7,5) * 0.5700448858423344 + O(7,7) * 0.4943528756111367);
    c_out[oi++] = (O(7,-7) * 0.2981060004427955 - O(7,-5) * 0.78125 -
                   O(7,-3) * 0.10364452469860624 +
                   O(7,-1) * 0.5385527481129402);
    c_out[oi++] = (O(7,1) * -0.4487939567607835 - O(7,3) * 0.5700448858423344 +
                   O(7,5) * 0.671875 - O(7,7) * 0.14905300022139775);
    c_out[oi++] = (O(7,-7) * -0.05846339666834283 +
                   O(7,-5) * 0.3983608994994363 -
                   O(7,-3) * 0.7927281808728639 + O(7,-1) * 0.4576818286211503);
    c_out[oi++] = (O(7,1) * -0.8562442974262661 + O(7,3) * 0.4943528756111367 -
                   O(7,5) * 0.14905300022139775 + O(7,7) * 0.015625);

    if (lmax < 8) {
      return;
    }
    c_out[oi++] = (O(8,1) * -0.8356088723200586 + O(8,3) * 0.516334738808072 -
                O(8,5) * 0.184877493221863 + O(8,7) * 0.03125);
    c_out[oi++] = (O(8,-7) * -0.109375 + O(8,-5) * 0.4621937330546575 - O(8,-3) * 0.774502108212108 +
                O(8,-1) * 0.4178044361600293);
    c_out[oi++] = (O(8,1) * -0.4576818286211503 - O(8,3) * 0.47134697278119864 +
                O(8,5) * 0.7088310138883598 - O(8,7) * 0.2567449488305466);
    c_out[oi++] = (O(8,-7) * 0.4621937330546575 - O(8,-5) * 0.703125 - O(8,-3) * 0.2181912506838897 +
                O(8,-1) * 0.4943528756111367);
    c_out[oi++] = (O(8,1) * -0.27421763710600383 - O(8,3) * 0.6051536478449089 -
                O(8,5) * 0.33802043207474897 + O(8,7) * 0.6665852814906732);
    c_out[oi++] = (O(8,-7) * -0.774502108212108 - O(8,-5) * 0.2181912506838897 +
                O(8,-3) * 0.265625 + O(8,-1) * 0.5310201708739509);
    c_out[oi++] = (O(8,1) * -0.1307281291459493 - O(8,3) * 0.38081430021731066 -
                O(8,5) * 0.5908647000371574 - O(8,7) * 0.6991205412874092);
    c_out[oi++] = (O(8,-7) * 0.4178044361600293 + O(8,-5) * 0.4943528756111367 +
                O(8,-3) * 0.5310201708739509 + O(8,-1) * 0.546875);
    c_out[oi++] = (O(8,0) * 0.2734375 + O(8,2) * 0.3921843874378479 + O(8,4) * 0.4113264556590057 +
                O(8,6) * 0.4576818286211503 + O(8,8) * 0.626706654240044);
    c_out[oi++] = (O(8,-8) * 0.8356088723200586 + O(8,-6) * 0.4576818286211503 +
                O(8,-4) * 0.27421763710600383 + O(8,-2) * 0.1307281291459493);
    c_out[oi++] = (O(8,0) * 0.3921843874378479 + O(8,2) * 0.5 + O(8,4) * 0.32775276505317236 -
                O(8,8) * 0.6991205412874092);
    c_out[oi++] = (O(8,-8) * -0.516334738808072 + O(8,-6) * 0.47134697278119864 +
                O(8,-4) * 0.6051536478449089 + O(8,-2) * 0.38081430021731066);
    c_out[oi++] = (O(8,0) * 0.4113264556590057 + O(8,2) * 0.32775276505317236 -
                O(8,4) * 0.28125 - O(8,6) * 0.7302075903467452 + O(8,8) * 0.3332926407453366);
    c_out[oi++] = (O(8,-8) * 0.184877493221863 - O(8,-6) * 0.7088310138883598 +
                O(8,-4) * 0.33802043207474897 + O(8,-2) * 0.5908647000371574);
    c_out[oi++] = (O(8,0) * 0.4576818286211503 - O(8,4) * 0.7302075903467452 + O(8,6) * 0.5 -
                O(8,8) * 0.0855816496101822);
    c_out[oi++] = (O(8,-8) * -0.03125 + O(8,-6) * 0.2567449488305466 - O(8,-4) * 0.6665852814906732 +
                O(8,-2) * 0.6991205412874092);
    c_out[oi++] = (O(8,0) * 0.626706654240044 - O(8,2) * 0.6991205412874092 +
                O(8,4) * 0.3332926407453366 - O(8,6) * 0.0855816496101822 + O(8,8) * 0.0078125);

    if (lmax < 9) {
      return;
    }
    c_out[oi++] = (O(9,0) * 0.6090493921755238 - O(9,2) * 0.6968469725305549 +
                O(9,4) * 0.3615761395439417 - O(9,6) * 0.11158481919598204 + O(9,8) * 0.016572815184059706);
    c_out[oi++] = (O(9,-8) * -0.0625 + O(9,-6) * 0.3156095293238149 - O(9,-4) * 0.6817945071647321 +
                O(9,-2) * 0.656993626300895);
    c_out[oi++] = (O(9,0) * 0.44314852502786806 - O(9,2) * 0.05633673867912483 - O(9,4) * 0.6723290616859425 +
                O(9,6) * 0.5683291712335379 - O(9,8) * 0.1594400908746762);
    c_out[oi++] = (O(9,-8) * 0.3156095293238149 - O(9,-6) * 0.71875 + O(9,-4) * 0.20252314682524564 +
                O(9,-2) * 0.5854685623498499);
    c_out[oi++] = (O(9,0) * 0.39636409043643195 + O(9,2) * 0.25194555463432966 - O(9,4) * 0.3921843874378479 -
                O(9,6) * 0.6051536478449089 + O(9,8) * 0.509312687906457);
    c_out[oi++] = (O(9,-8) * -0.6817945071647321 + O(9,-6) * 0.20252314682524564 + O(9,-4) * 0.5625 +
                O(9,-2) * 0.4215855488510013);
    c_out[oi++] = (O(9,0) * 0.3754879637718099 + O(9,2) * 0.42961647140211 + O(9,4) * 0.13799626353637262 -
                O(9,6) * 0.2981060004427955 - O(9,8) * 0.7526807559068452);
    c_out[oi++] = (O(9,-8) * 0.656993626300895 + O(9,-6) * 0.5854685623498499 + O(9,-4) * 0.4215855488510013 +
                O(9,-2) * 0.21875);
    c_out[oi++] = (O(9,0) * 0.36685490255855924 + O(9,2) * 0.5130142237306876 + O(9,4) * 0.4943528756111367 +
                O(9,6) * 0.4576818286211503 + O(9,8) * 0.38519665736315783);
    c_out[oi++] = (O(9,-9) * -0.6090493921755238 - O(9,-7) * 0.44314852502786806 - O(9,-5) * 0.39636409043643195 -
                O(9,-3) * 0.3754879637718099 - O(9,-1) * 0.36685490255855924);
    c_out[oi++] = (O(9,1) * 0.0546875 + O(9,3) * 0.16792332234534904 + O(9,5) * 0.2954323500185787 +
                O(9,7) * 0.4624247721758373 + O(9,9) * 0.8171255055356398);
    c_out[oi++] = (O(9,-9) * 0.6968469725305549 + O(9,-7) * 0.05633673867912483 - O(9,-5) * 0.25194555463432966 -
                O(9,-3) * 0.42961647140211 - O(9,-1) * 0.5130142237306876);
    c_out[oi++] = (O(9,1) * 0.16792332234534904 + O(9,3) * 0.453125 + O(9,5) * 0.577279787559724 +
                O(9,7) * 0.387251054106054 - O(9,9) * 0.5322256665703469);
    c_out[oi++] = (O(9,-9) * -0.3615761395439417 + O(9,-7) * 0.6723290616859425 + O(9,-5) * 0.3921843874378479 -
                O(9,-3) * 0.13799626353637262 - O(9,-1) * 0.4943528756111367);
    c_out[oi++] = (O(9,1) * 0.2954323500185787 + O(9,3) * 0.577279787559724 + O(9,5) * 0.140625 -
                O(9,7) * 0.7162405240429014 + O(9,9) * 0.21608307321780204);
    c_out[oi++] = (O(9,-9) * 0.11158481919598204 - O(9,-7) * 0.5683291712335379 + O(9,-5) * 0.6051536478449089 +
                O(9,-3) * 0.2981060004427955 - O(9,-1) * 0.4576818286211503);
    c_out[oi++] = (O(9,1) * 0.4624247721758373 + O(9,3) * 0.387251054106054 - O(9,5) * 0.7162405240429014 +
                O(9,7) * 0.34765625 - O(9,9) * 0.048317644050206957);
    c_out[oi++] = (O(9,-9) * -0.016572815184059706 +
                   O(9,-7) * 0.1594400908746762 - O(9,-5) * 0.509312687906457 +
                   O(9,-3) * 0.7526807559068452 -
                   O(9,-1) * 0.38519665736315783);
    c_out[oi++] = (O(9,1) * 0.8171255055356398 - O(9,3) * 0.5322256665703469 +
                   O(9,5) * 0.21608307321780204 -
                   O(9,7) * 0.048317644050206957 + O(9,9) * 0.00390625);
  }

  static void ConvolveCosTheta(int lmax, List<Spectrum> c_in,
                               List<Spectrum> c_out) {
    const List<double> c_costheta = const [ 0.8862268925, 1.0233267546,
            0.4954159260, 0.0000000000, -0.1107783690, 0.0000000000,
            0.0499271341, 0.0000000000, -0.0285469331, 0.0000000000,
            0.0185080823, 0.0000000000, -0.0129818395, 0.0000000000,
            0.0096125342, 0.0000000000, -0.0074057109, 0.0000000000 ];
    for (int l = 0; l <= lmax; ++l) {
      for (int m = -l; m <= l; ++m) {
        int o = Index(l, m);
        if (l < 18) {
          c_out[o] = c_in[o] * (_lambda(l) * c_costheta[l]);
        } else {
          c_out[o].set(0.0);
        }
      }
    }
  }

  static void ConvolvePhong(int lmax, double n, List<Spectrum> c_in,
                            List<Spectrum> c_out) {
    for (int l = 0; l <= lmax; ++l) {
      double c_phong = Math.exp(-(l * l) / (2.0 * n));
      for (int m = -l; m <= l; ++m) {
        int o = Index(l, m);
        c_out[o] = c_in[o] * (_lambda(l) * c_phong);
      }
    }
  }

  static void ComputeDiffuseTransfer(Point p, Normal n, double rayEpsilon,
                                     Scene scene, RNG rng, int nSamples,
                                     int lmax, List<Spectrum> c_transfer) {
    List<int> scramble = [ rng.randomUint(), rng.randomUint() ];
    List<double> Ylm = new List<double>(Terms(lmax));
    List<double> u = [0.0, 0.0];
    for (int i = 0; i < nSamples; ++i) {
      // Sample _i_th direction and compute estimate for transfer coefficients
      Sample02(i, scramble, u);
      Vector w = UniformSampleSphere(u[0], u[1]);
      double pdf = UniformSpherePdf();
      if (Vector.Dot(w, n) > 0.0 && !scene.intersectP(new Ray(p, w, rayEpsilon))) {
        // Accumulate contribution of direction $\w{}$ to transfer coefficients
        Evaluate(w, lmax, Ylm);
        for (int j = 0, len = Terms(lmax); j < len; ++j) {
          c_transfer[j] += new Spectrum(Ylm[j] * Vector.AbsDot(w, n)) /
                          (pdf * nSamples);
        }
      }
    }
  }

  static void ComputeTransferMatrix(Point p, double rayEpsilon, Scene scene,
                                    RNG rng, int nSamples, int lmax,
                                    List<Spectrum> T) {
    for (int i = 0, len = Terms(lmax) * Terms(lmax); i < len; ++i) {
      T[i] = new Spectrum(0.0);
    }
    List<int> scramble = [ rng.randomUint(), rng.randomUint() ];
    List<double> Ylm = new List<double>(Terms(lmax));
    List<double> u = [0.0, 0.0];
    for (int i = 0; i < nSamples; ++i) {
      // Compute Monte Carlo estimate of $i$th sample for transfer matrix
      Sample02(i, scramble, u);
      Vector w = UniformSampleSphere(u[0], u[1]);
      double pdf = UniformSpherePdf();
      if (!scene.intersectP(new Ray(p, w, rayEpsilon))) {
        // Update transfer matrix for unoccluded direction
        Evaluate(w, lmax, Ylm);
        for (int j = 0, nj = Terms(lmax); j < nj; ++j) {
          for (int k = 0, nk = Terms(lmax); k < nk; ++k) {
            T[j * Terms(lmax) + k] += new Spectrum((Ylm[j] * Ylm[k]) /
                                                     (pdf * nSamples));
          }
        }
      }
    }
  }

  static void ComputeBSDFMatrix(Spectrum Kd, Spectrum Ks, double roughness,
                                RNG rng, int nSamples, int lmax,
                                List<Spectrum> B) {
    for (int i = 0; i < Terms(lmax) * Terms(lmax); ++i) {
      B[i] = new Spectrum(0.0);
    }

    // Create _BSDF_ for computing BSDF transfer matrix
    DifferentialGeometry dg = new DifferentialGeometry().set(
                                                     new Point(0.0, 0.0, 0.0),
                                                     new Vector(1.0, 0.0, 0.0),
                                                     new Vector(0.0, 1.0, 0.0),
                                                     new Normal(0.0, 0.0, 0.0),
                                                     new Normal(0.0, 0.0, 0.0),
                                                     0.0, 0.0, null);

    BSDF bsdf = new BSDF(dg, new Normal(0.0, 0.0, 1.0));

    bsdf.add(new Lambertian(new Spectrum.from(Kd)));

    Fresnel fresnel = new FresnelDielectric(1.5, 1.0);

    bsdf.add(new Microfacet(Ks, fresnel, new Blinn(1.0 / roughness)));

    // Precompute directions w{} and SH values for directions
    Float32List Ylm = new Float32List(Terms(lmax) * nSamples);
    List<Vector> w = new List<Vector>(nSamples);

    List<int> scramble = [ rng.randomUint(), rng.randomUint() ];
    List<double> u = [0.0, 0.0];
    for (int i = 0; i < nSamples; ++i) {
      Sample02(i, scramble, u);
      w[i] = UniformSampleSphere(u[0], u[1]);
      Evaluate(w[i], lmax, Ylm, Terms(lmax) * i);
    }

    // Compute double spherical integral for BSDF matrix
    for (int osamp = 0; osamp < nSamples; ++osamp) {
      Vector wo = w[osamp];
      for (int isamp = 0; isamp < nSamples; ++isamp) {
        Vector wi = w[isamp];
        // Update BSDF matrix elements for sampled directions
        Spectrum f = bsdf.f(wo, wi);
        if (!f.isBlack()) {
          double pdf = UniformSpherePdf() * UniformSpherePdf();
          f *= (Vector.CosTheta(wi)).abs() / (pdf * nSamples * nSamples);
          for (int i = 0; i < Terms(lmax); ++i) {
            for (int j = 0; j < Terms(lmax); ++j) {
                B[i * Terms(lmax) + j].add(f * (Ylm[isamp * Terms(lmax) + j] *
                                             Ylm[osamp * Terms(lmax) + i]));
            }
          }
        }
      }
    }
  }

  static void MatrixVectorMultiply(List<Spectrum> M, List<Spectrum> v,
                                   List<Spectrum> vout, int lmax) {
    for (int i = 0, len = Terms(lmax); i < len; ++i) {
      vout[i] = new Spectrum(0.0);
      for (int j = 0; j < Terms(lmax); ++j) {
        vout[i] += v[j] * M[Terms(lmax) * i + j];
      }
    }
  }

  static void _legendrep(double x, int lmax, List<double> out,
                         [int outIndex = 0]) {
    P(int l, int m) {
      int index = (outIndex + Index(l, m)).toInt();
      return out[index];
    }

    // Compute m=0 Legendre values using recurrence
    out[Index(0, 0)] = 1.0;
    out[Index(1, 0)] = x;
    for (int l = 2; l <= lmax; ++l) {
      out[outIndex + Index(l, 0)] = ((2.0 * l - 1.0) * x * P(l - 1, 0) -
                            (l - 1) * P(l - 2, 0)) / l;
      assert(!P(l, 0).isNaN);
      assert(P(l, 0).isFinite);
    }

    // Compute $m=l$ edge using Legendre recurrence
    double neg = -1.0;
    double dfact = 1.0;
    double xroot = Math.sqrt(Math.max(0.0, 1.0 - x * x));
    double xpow = xroot;
    for (int l = 1; l <= lmax; ++l) {
      out[outIndex + Index(l, l)] = neg * dfact * xpow;
      assert(!P(l, l).isNaN);
      assert(P(l, l).isFinite);
      neg *= -1.0;      // neg = (-1)^l
      dfact *= 2 * l + 1; // dfact = (2*l-1)!!
      xpow *= xroot;    // xpow = powf(1.f - x*x, double(l) * 0.5f);
    }

    // Compute $m=l-1$ edge using Legendre recurrence
    for (int l = 2; l <= lmax; ++l) {
      out[outIndex + Index(l, l - 1)] = x * (2.0 * l - 1.0) * P(l - 1, l - 1);
      assert(!P(l, l - 1).isNaN);
      assert(P(l, l - 1).isFinite);
    }

    // Compute $m=1, \ldots, l-2$ values using Legendre recurrence
    for (int l = 3; l <= lmax; ++l) {
      for (int m = 1; m <= l - 2; ++m) {
        out[outIndex + Index(l, m)] =
            ((2.0 * (l - 1.0) + 1.0) * x * P(l - 1, m) -
             (l - 1.0 + m) * P(l - 2, m)) / (l - m);
        assert(!P(l, m).isNaN);
        assert(P(l, m).isFinite);
      }
    }
  }

  static double _K(int l, int m) {
    return Math.sqrt((2.0 * l + 1.0) * INV_FOURPI * _divfact(l, m));
  }

  static double _divfact(int a, int b) {
    if (b == 0) {
      return 1.0;
    }
    double fa = a.toDouble();
    double fb = b.toDouble().abs();
    double v = 1.0;
    for (double x = fa - fb + 1.0; x <= fa + fb; x += 1.0) {
      v *= x;
    }
    return 1.0 / v;
  }

  /**
   * n!! = 1 if n==0 or 1, otherwise n * (n-2)!!
   */
  static double _dfact(double v) {
    if (v <= 1.0) {
      return 1.0;
    }
    return v * _dfact(v - 2.0);
  }

  static double _fact(double v) {
    if (v <= 1.0) {
      return 1.0;
    }
    return v * _fact(v - 1.0);
  }

  static void _sinCosIndexed(double s, double c, int n,
                             List<double> sout, List<double> cout) {
    double si = 0.0;
    double ci = 1.0;
    for (int i = 0; i < n; ++i) {
      // Compute $\sin{}i\phi$ and $\cos{}i\phi$ using recurrence
      sout[i] = si;
      cout[i] = ci;
      double oldsi = si;
      si = si * c + ci * s;
      ci = ci * c - oldsi * s;
    }
  }

  static void _toZYZ(Matrix4x4 m, List<double> alpha, List<double> beta,
                     List<double> gamma) {
    double sy = Math.sqrt(m[9] * m[9] + m[8] * m[8]);
    if (sy > 16 * FLT_EPSILON) {
      gamma[0] = -Math.atan2(m[6], -m[2]);
      beta[0]  = -Math.atan2(sy, m[10]);
      alpha[0] = -Math.atan2(m[9], m[8]);
    } else {
      gamma[0] =  0.0;
      beta[0]  = -Math.atan2(sy, m[10]);
      alpha[0] = -Math.atan2(-m[4], m[5]);
    }
  }


  static double _lambda(num l) =>
      Math.sqrt((4.0 * Math.PI) / (2.0 * l + 1.0));
}
