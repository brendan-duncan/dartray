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
part of surface_integrators;

class UseProbesIntegrator extends SurfaceIntegrator {
  UseProbesIntegrator(String filename) {
    Completer completer = new Completer();
    if (filename.isNotEmpty) {
      // Read precomputed radiance probe values from file
      ResourceManager.RequestFile(filename, completer.future).then((bytes) {

        List<double> fp;
        if (bytes is List<int>) {
          fp = ReadFloatFile(bytes, filename);
          ResourceManager.WriteFile(filename, fp);
        } else if (bytes is List<double>) {
          fp = bytes;
        }

        int fi = 0;

        lmax = fp[fi++].toInt();
        includeDirectInProbes = fp[fi++].toInt();
        includeIndirectInProbes = fp[fi++].toInt();

        nProbes[0] = fp[fi++].toInt();
        nProbes[1] = fp[fi++].toInt();
        nProbes[2] = fp[fi++].toInt();

        bbox.pMin.x = fp[fi++];
        bbox.pMin.y = fp[fi++];
        bbox.pMin.z = fp[fi++];
        bbox.pMax.x = fp[fi++];
        bbox.pMax.y = fp[fi++];
        bbox.pMax.z = fp[fi++];

        final int lmax_terms = SphericalHarmonics.Terms(lmax);
        final int numProbes = nProbes[0] * nProbes[1] * nProbes[2];

        c_in = new List<Spectrum>(lmax_terms * numProbes);

        int offset = 0;
        for (int i = 0; i < numProbes; ++i) {
          for (int j = 0; j < lmax_terms; ++j) {
            c_in[offset++] = new Spectrum()..setList(fp, fi);
            fi += Spectrum.NumSamples();
          }
        }

        completer.complete();
      });
    }
  }

  void requestSamples(Sampler sampler, Sample sample, Scene scene) {
    // Allocate and request samples for sampling all lights
    int nLights = scene.lights.length;
    lightSampleOffsets = new List<LightSampleOffsets>(nLights);
    bsdfSampleOffsets = new List<BSDFSampleOffsets>(nLights);
    for (int i = 0; i < nLights; ++i) {
      Light light = scene.lights[i];
      int nSamples = light.nSamples;
      if (sampler != null) {
        nSamples = sampler.roundSize(nSamples);
      }
      lightSampleOffsets[i] = new LightSampleOffsets(nSamples, sample);
      bsdfSampleOffsets[i] = new BSDFSampleOffsets(nSamples, sample);
    }
  }

  Spectrum Li(Scene scene, Renderer renderer,
              RayDifferential ray, Intersection isect,
              Sample sample, RNG rng) {
    Spectrum L = new Spectrum(0.0);
    Vector wo = -ray.direction;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Evaluate BSDF at hit point
    BSDF bsdf = isect.getBSDF(ray);
    Point p = bsdf.dgShading.p;
    Normal n = bsdf.dgShading.nn;

    // Compute reflection for radiance probes integrator
    if (includeDirectInProbes == 0) {
      L += Integrator.UniformSampleAllLights(scene, renderer, p, n,
                                             wo, isect.rayEpsilon, ray.time,
                                             bsdf, sample, rng,
                                             lightSampleOffsets,
                                             bsdfSampleOffsets);
    }

    // Compute reflected lighting using radiance probes

    // Compute probe coordinates and offsets for lookup point
    Vector offset = bbox.offset(p);
    double voxx = (offset.x * nProbes[0]) - 0.5;
    double voxy = (offset.y * nProbes[1]) - 0.5;
    double voxz = (offset.z * nProbes[2]) - 0.5;
    int vx = voxx.floor();
    int vy = voxy.floor();
    int vz = voxz.floor();
    double dx = voxx - vx;
    double dy = voxy - vy;
    double dz = voxz - vz;

    // Get radiance probe coefficients around lookup point
    int b000 = c_inXYZ(lmax, vx,   vy,   vz);
    int b100 = c_inXYZ(lmax, vx+1, vy,   vz);
    int b010 = c_inXYZ(lmax, vx,   vy+1, vz);
    int b110 = c_inXYZ(lmax, vx+1, vy+1, vz);
    int b001 = c_inXYZ(lmax, vx,   vy,   vz+1);
    int b101 = c_inXYZ(lmax, vx+1, vy,   vz+1);
    int b011 = c_inXYZ(lmax, vx,   vy+1, vz+1);
    int b111 = c_inXYZ(lmax, vx+1, vy+1, vz+1);

    // Compute incident radiance from radiance probe coefficients
    int lmax_terms = SphericalHarmonics.Terms(lmax);
    List<Spectrum> c_inp = new List<Spectrum>(lmax_terms);
    for (int i = 0; i < lmax_terms; ++i) {
      // Do trilinear interpolation to compute SH coefficients at point
      Spectrum c00 = Lerp(dx, c_in[b000 + i], c_in[b100 + i]);
      Spectrum c10 = Lerp(dx, c_in[b010 + i], c_in[b110 + i]);
      Spectrum c01 = Lerp(dx, c_in[b001 + i], c_in[b101 + i]);
      Spectrum c11 = Lerp(dx, c_in[b011 + i], c_in[b111 + i]);
      Spectrum c0 = Lerp(dy, c00, c10);
      Spectrum c1 = Lerp(dy, c01, c11);
      c_inp[i] = Lerp(dz, c0, c1);
    }

    // Convolve incident radiance to compute irradiance function
    List<Spectrum> c_E = new List<Spectrum>(lmax_terms);
    SphericalHarmonics.ConvolveCosTheta(lmax, c_inp, c_E);

    // Evaluate irradiance function and accumulate reflection
    Spectrum rho = bsdf.rho2(wo, rng, BSDF_ALL_REFLECTION);
    Float32List Ylm = new Float32List(lmax_terms);
    SphericalHarmonics.Evaluate(Vector.FaceForward(n, wo), lmax, Ylm);
    Spectrum E = new Spectrum(0.0);
    for (int i = 0; i < lmax_terms; ++i) {
      E += c_E[i] * Ylm[i];
    }
    L += rho * INV_PI * E.clamp();

    return L;
  }

  int c_inXYZ(int lmax, int vx, int vy, int vz) {
    vx = vx.clamp(0, nProbes[0] - 1);
    vy = vy.clamp(0, nProbes[1] - 1);
    vz = vz.clamp(0, nProbes[2] - 1);
    int offset = vx + vy * nProbes[0] + vz * nProbes[0] * nProbes[1];
    return SphericalHarmonics.Terms(lmax) * offset;
  }

  static UseProbesIntegrator Create(ParamSet paramSet) {
    String filename = paramSet.findOneFilename("filename", "probes.out");
    return new UseProbesIntegrator(filename);
  }

  BBox bbox = new BBox();
  int lmax;
  int includeDirectInProbes;
  int includeIndirectInProbes;
  List<int> nProbes = [0, 0, 0];
  List<Spectrum> c_in;
  List<LightSampleOffsets> lightSampleOffsets;
  List<BSDFSampleOffsets> bsdfSampleOffsets;
}
