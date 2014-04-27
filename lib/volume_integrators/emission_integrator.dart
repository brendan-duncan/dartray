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
part of volume_integrators;

class EmissionIntegrator extends VolumeIntegrator {
  EmissionIntegrator(this.stepSize);

  static EmissionIntegrator Create(ParamSet params) {
    double stepSize = params.findOneFloat('stepsize', 1.0);
    return new EmissionIntegrator(stepSize);
  }

  void requestSamples(Sampler sampler, Sample sample, Scene scene) {
    tauSampleOffset = sample.add1D(1);
    scatterSampleOffset = sample.add1D(1);
  }

  Spectrum Li(Scene scene, Renderer renderer,
              RayDifferential ray, Sample sample, RNG rng, Spectrum T) {
    VolumeRegion vr = scene.volumeRegion;
    //Assert(sample != NULL);
    List<double> t0 = [0.0];
    List<double> t1 = [0.0];

    if (vr == null || !vr.intersectP(ray, t0, t1) || (t1[0] - t0[0]) == 0.0) {
      T.set(1.0);
      return new Spectrum(0.0);
    }

    // Do emission-only volume integration in _vr_
    Spectrum Lv = new Spectrum(0.0);

    // Prepare for volume integration stepping
    int nSamples = ((t1[0] - t0[0]) / stepSize).ceil();
    double step = (t1[0] - t0[0]) / nSamples;
    Spectrum Tr = new Spectrum(1.0);
    Point p = ray.pointAt(t0[0]);
    Point pPrev = new Point();
    Vector w = -ray.direction;

    t0[0] += sample.oneD[scatterSampleOffset][0] * step;

    for (int i = 0; i < nSamples; ++i, t0[0] += step) {
      // Advance to sample at _t0_ and update _T_
      pPrev = p;
      p = ray.pointAt(t0[0]);
      Ray tauRay = new Ray(pPrev, p - pPrev, 0.0, 1.0, ray.time, ray.depth);
      Spectrum stepTau = vr.tau(tauRay, 0.5 * stepSize, rng.randomFloat());
      Tr *= (-stepTau).exp();

      // Possibly terminate ray marching if transmittance is small
      if (Tr.luminance() < 1.0e-3) {
        const double continueProb = 0.5;
        if (rng.randomFloat() > continueProb) {
          Tr.set(0.0);
          break;
        }
        Tr /= continueProb;
      }

      // Compute emission-only source term at _p_
      Lv += Tr * vr.Lve(p, w, ray.time);
    }

    T.copy(Tr);
    return Lv * step;
  }

  Spectrum transmittance(Scene scene, Renderer renderer,
                         RayDifferential ray, Sample sample, RNG rng) {
    if (scene.volumeRegion == null) {
      return new Spectrum(1.0);
    }

    double step;
    double offset;
    if (sample != null) {
      step = stepSize;
      offset = sample.oneD[tauSampleOffset][0];
    } else {
      step = 4.0 * stepSize;
      offset = rng.randomFloat();
    }

    Spectrum tau = scene.volumeRegion.tau(ray, step, offset);

    return (-tau).exp();
  }

  double stepSize;
  int tauSampleOffset;
  int scatterSampleOffset;
}
