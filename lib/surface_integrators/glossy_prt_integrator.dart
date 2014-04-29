/****************************************************************************
 * Copyright (C) 2014 by Brendan Duncan.                                    *
 *                                                                          *
 * This file is part of DartRay.                                            *
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the 'License');          *
 * you may not use this file except in compliance with the License.         *
 * You may obtain a copy of the License at                                  *
 *                                                                          *
 * http://www.apache.org/licenses/LICENSE-2.0                               *
 *                                                                          *
 * Unless required by applicable law or agreed to in writing, software      *
 * distributed under the License is distributed on an 'AS IS' BASIS,        *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 * See the License for the specific language governing permissions and      *
 * limitations under the License.                                           *
 *                                                                          *
 * This project is based on PBRT v2 ; see http://www.pbrt.org               *
 * pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.  *
 ****************************************************************************/
part of surface_integrators;

class GlossyPRTIntegrator extends SurfaceIntegrator {
  GlossyPRTIntegrator(this.Kd, this.Ks, this.roughness, this.lmax, int ns)
      : nSamples = RoundUpPow2(ns);

  void preprocess(Scene scene, Camera camera, Renderer renderer) {
    // Project direct lighting into SH for GlossyPRTIntegrator
    Stopwatch timer = new Stopwatch()..start();
    LogDebug('STARTING GlossyPRT Preprocess');
    BBox bbox = scene.worldBound;
    Point p = bbox.center;
    RNG rng = new RNG();

    int lmax_terms = SphericalHarmonics.Terms(lmax);

    c_in = new List<Spectrum>(lmax_terms);
    for (int i = 0, len = c_in.length; i < len; ++i) {
      c_in[i] = new Spectrum(0.0);
    }

    SphericalHarmonics.ProjectIncidentDirectRadiance(p, 0.0, camera.shutterOpen,
                                                     scene, false, lmax, rng,
                                                     c_in);

    // Compute glossy BSDF matrix for PRT
    B = new List<Spectrum>(lmax_terms * lmax_terms);
    SphericalHarmonics.ComputeBSDFMatrix(Kd, Ks, roughness, rng, 1024, lmax, B);

    LogDebug('FINISH GlossyPRT Preprocess: ${timer.elapsed}');
  }

  void requestSamples(Sampler sampler, Sample sample, Scene scene) {
  }

  Spectrum Li(Scene scene, Renderer renderer, RayDifferential ray,
              Intersection isect, Sample sample, RNG rng) {
    Spectrum L = new Spectrum(0.0);
    Vector wo = -ray.direction;

    final int lmax_terms = SphericalHarmonics.Terms(lmax);

    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Evaluate BSDF at hit point
    BSDF bsdf = isect.getBSDF(ray);
    Point p = bsdf.dgShading.p;

    // Compute reflected radiance with glossy PRT at point

    // Compute SH radiance transfer matrix at point and SH coefficients
    List<Spectrum> c_t = new List<Spectrum>(lmax_terms);
    List<Spectrum> T = new List<Spectrum>(lmax_terms * lmax_terms);

    SphericalHarmonics.ComputeTransferMatrix(p, isect.rayEpsilon, scene, rng,
                                             nSamples, lmax, T);

    SphericalHarmonics.MatrixVectorMultiply(T, c_in, c_t, lmax);

    // Rotate incident SH lighting to local coordinate frame
    Vector r1 = bsdf.localToWorld(new Vector(1.0, 0.0, 0.0));
    Vector r2 = bsdf.localToWorld(new Vector(0.0, 1.0, 0.0));
    Normal nl = new Normal.from(bsdf.localToWorld(new Vector(0.0, 0.0, 1.0)));
    Matrix4x4 rot = new Matrix4x4.values(
                            r1.x, r2.x, nl.x, 0.0,
                            r1.y, r2.y, nl.y, 0.0,
                            r1.z, r2.z, nl.z, 0.0,
                            0.0,  0.0,  0.0,  1.0);

    List<Spectrum> c_l = new List<Spectrum>(lmax_terms);
    for (int i = 0, len = c_l.length; i < len; ++i) {
      c_l[i] = new Spectrum(0.0);
    }
    SphericalHarmonics.Rotate(c_t, c_l, rot, lmax);

    // Compute final coefficients c_o using BSDF matrix
    List<Spectrum> c_o = new List<Spectrum>(lmax_terms);
    SphericalHarmonics.MatrixVectorMultiply(B, c_l, c_o, lmax);

    // Evaluate outgoing radiance function for wo and add to L
    Vector woLocal = bsdf.worldToLocal(wo);

    Float32List Ylm = new Float32List(lmax_terms);
    SphericalHarmonics.Evaluate(woLocal, lmax, Ylm);

    Spectrum Li = new Spectrum(0.0);

    for (int i = 0; i < lmax_terms; ++i) {
      Li += c_o[i] * Ylm[i];
    }

    L += Li.clamp();

    return L;
  }

  static GlossyPRTIntegrator Create(ParamSet params) {
    int lmax = params.findOneInt('lmax', 4);
    int ns = params.findOneInt('nsamples', 4096);
    Spectrum Kd = params.findOneSpectrum('Kd', new Spectrum(0.5));
    Spectrum Ks = params.findOneSpectrum('Ks', new Spectrum(0.25));
    double roughness = params.findOneFloat('roughness', 0.1);
    return new GlossyPRTIntegrator(Kd, Ks, roughness, lmax, ns);
  }

  Spectrum Kd;
  Spectrum Ks;
  double roughness;
  int lmax;
  int nSamples;
  List<Spectrum> c_in;
  List<Spectrum> B;
}
