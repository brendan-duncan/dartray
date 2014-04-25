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

class WhittedIntegrator extends SurfaceIntegrator {
  WhittedIntegrator(this.maxDepth);

  Spectrum Li(Scene scene, Renderer renderer,
                 RayDifferential ray, Intersection isect, Sample sample,
                 RNG rng) {
    Spectrum L = new Spectrum(0.0);
    // Compute emitted and reflected light at ray intersection point

    // Evaluate BSDF at hit point
    BSDF bsdf = isect.getBSDF(ray);

    // Initialize common variables for Whitted integrator
    Point p = bsdf.dgShading.p;
    Normal n = bsdf.dgShading.nn;
    Vector wo = -ray.direction;

    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Add contribution of each light source
    for (int i = 0; i < scene.lights.length; ++i) {
      Vector wi = new Vector();
      List<double> pdf = [0.0];
      VisibilityTester visibility = new VisibilityTester();
      Spectrum Li = scene.lights[i].sampleLAtPoint(p, isect.rayEpsilon,
          new LightSample.random(rng), ray.time, wi, pdf, visibility);

      if (Li.isBlack() || pdf[0] == 0.0) {
        continue;
      }

      Spectrum f = bsdf.f(wo, wi);

      if (!f.isBlack() && visibility.unoccluded(scene)) {
          L += f * Li * Vector.AbsDot(wi, n) *
               visibility.transmittance(scene, renderer, sample, rng) / pdf[0];
      }
    }

    if (ray.depth + 1 < maxDepth) {
      // Trace rays for specular reflection and refraction
      L += Integrator.SpecularReflect(ray, bsdf, rng, isect, renderer, scene,
                                      sample);
      L += Integrator.SpecularTransmit(ray, bsdf, rng, isect, renderer, scene,
                                       sample);
    }

    return L;
  }

  static WhittedIntegrator Create(ParamSet params) {
    int maxDepth = params.findOneInt("maxdepth", 5);
    return new WhittedIntegrator(maxDepth);
  }

  int maxDepth;
}
