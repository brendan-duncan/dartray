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

class DipoleSubsurfaceIntegrator extends SurfaceIntegrator {
  DipoleSubsurfaceIntegrator(int mdepth, double merror, double mindist,
                             String fn) {
    maxSpecularDepth = mdepth;
    maxError = merror;
    minSampleDist = mindist;
    filename = fn;
    octree = null;

    // Initialize _SurfacePoint_s from file
    Completer completer = new Completer();
    if (filename.isNotEmpty) {
      ResourceManager.RequestFile(filename, completer.future).then((bytes) {
        if (bytes is List<int>) {
          List<double> floats = ReadFloatFile(bytes, filename);
          ResourceManager.WriteFile(filename, floats);
        }
        completer.complete();
      });
    }
  }

  Spectrum Li(Scene scene, Renderer renderer,
              RayDifferential ray, Intersection isect, Sample sample,
              RNG rng) {
    Spectrum L = new Spectrum(0.0);
    Vector wo = -ray.direction;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Evaluate BSDF at hit point
    BSDF bsdf = isect.getBSDF(ray);
    Point p = bsdf.dgShading.p;
    Normal n = bsdf.dgShading.nn;

    // Evaluate BSSRDF and possibly compute subsurface scattering
    BSSRDF bssrdf = isect.getBSSRDF(ray);

    if (bssrdf != null && octree != null) {
      Spectrum sigma_a  = bssrdf.sigma_a;
      Spectrum sigmap_s = bssrdf.sigma_prime_s;
      Spectrum sigmap_t = sigmap_s + sigma_a;

      if (!sigmap_t.isBlack()) {
        // Use hierarchical integration to evaluate reflection from dipole model
        Stats.SUBSURFACE_STARTED_OCTREE_LOOKUP(p);
        _DiffusionReflectance Rd = new _DiffusionReflectance(sigma_a, sigmap_s,
                                                             bssrdf.eta);
        Spectrum Mo = octree.Mo(octreeBounds, p, Rd, maxError);

        FresnelDielectric fresnel = new FresnelDielectric(1.0, bssrdf.eta);

        Spectrum Ft = Spectrum.ONE - fresnel.evaluate(Vector.AbsDot(wo, n));
        double Fdt = 1.0 - Fdr(bssrdf.eta);

        L += (Ft * INV_PI) * (Mo * Fdt);

        Stats.SUBSURFACE_FINISHED_OCTREE_LOOKUP();
      }
    }

    L += Integrator.UniformSampleAllLights(scene, renderer, p, n,
                                           wo, isect.rayEpsilon, ray.time,
                                           bsdf, sample, rng,
                                           lightSampleOffsets,
                                           bsdfSampleOffsets);

    if (ray.depth < maxSpecularDepth) {
      // Trace rays for specular reflection and refraction
      L += Integrator.SpecularReflect(ray, bsdf, rng, isect, renderer, scene,
                                      sample);

      L += Integrator.SpecularTransmit(ray, bsdf, rng, isect, renderer, scene,
                                       sample);
    }

    return L;
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

  void preprocess(Scene scene, Camera camera, Renderer renderer) {
    if (scene.lights.isEmpty) {
      return;
    }

    Stopwatch timer = new Stopwatch()..start();
    LogInfo('STARTING DipoleSubsurface Preprocessing. This may take a while.');

    List<SurfacePoint> pts = [];

    // Get _SurfacePoint_s for translucent objects in scene
    if (filename.isNotEmpty) {
      List<double> fpts = ResourceManager.GetResource(filename);
      if (fpts != null) {
        if ((fpts.length % 8) != 0) {
          LogError("Excess values (${fpts.length}) in points "
                   "file \"${filename}\"");
        }

        for (int i = 0; i < fpts.length; i += 8) {
          pts.add(new SurfacePoint(new Point(fpts[i], fpts[i + 1], fpts[i + 2]),
                              new Normal(fpts[i + 3], fpts[i + 4], fpts[i + 5]),
                              fpts[i + 6], fpts[i + 7]));
        }
      }
    }

    if (pts.isEmpty) {
      Point pCamera = camera.cameraToWorld.transformPoint(camera.shutterOpen,
                                                          Point.ZERO);

      SurfacePointsRenderer.FindPoissonPointDistribution(pCamera,
                                                        camera.shutterOpen,
                                                        minSampleDist, scene,
                                                        pts);
    }

    // Compute irradiance values at sample points
    RNG rng = new RNG();
    Stats.SUBSURFACE_STARTED_COMPUTING_IRRADIANCE_VALUES();

    List<double> lpos = [0.0, 0.0];

    for (int i = 0; i < pts.length; ++i) {
      SurfacePoint sp = pts[i];
      Spectrum E = new Spectrum(0.0);

      for (int j = 0; j < scene.lights.length; ++j) {
        // Add irradiance from light at point
        Light light = scene.lights[j];
        Spectrum Elight = new Spectrum(0.0);
        int nSamples = RoundUpPow2(light.nSamples);
        List<int> scramble = [rng.randomUint(), rng.randomUint()];
        int compScramble = rng.randomUint();
        for (int s = 0; s < nSamples; ++s) {
          Sample02(s, scramble, lpos);
          double lcomp = VanDerCorput(s, compScramble);
          LightSample ls = new LightSample(lpos[0], lpos[1], lcomp);
          Vector wi = new Vector();
          List<double> lightPdf = [0.0];
          VisibilityTester visibility = new VisibilityTester();
          Spectrum Li = light.sampleLAtPoint(sp.p, sp.rayEpsilon, ls,
                                             camera.shutterOpen, wi, lightPdf,
                                             visibility);
          if (Vector.Dot(wi, sp.n) <= 0.0) {
            continue;
          }

          if (Li.isBlack() || lightPdf[0] == 0.0) {
            continue;
          }

          Li *= visibility.transmittance(scene, renderer, null, rng);
          if (visibility.unoccluded(scene)) {
            Elight += (Li * Vector.AbsDot(wi, sp.n) / lightPdf[0]);
          }
        }

        E += Elight / nSamples;
      }

      irradiancePoints.add(new _IrradiancePoint(sp, E));

      Stats.SUBSURFACE_COMPUTED_IRRADIANCE_AT_POINT(sp, E);
    }

    Stats.SUBSURFACE_FINISHED_COMPUTING_IRRADIANCE_VALUES();

    // Create octree of clustered irradiance samples
    octree = new _SubsurfaceOctreeNode();

    octreeBounds = new BBox();
    for (int i = 0; i < irradiancePoints.length; ++i) {
      octreeBounds = BBox.UnionPoint(octreeBounds, irradiancePoints[i].p);
    }

    for (int i = 0; i < irradiancePoints.length; ++i) {
      octree.insert(octreeBounds, irradiancePoints[i]);
    }

    octree.initHierarchy();
    LogInfo('FINISHED DipoleSubsurface Preprocessing: ${timer.elapsed}');
  }

  int maxSpecularDepth;
  double maxError;
  double minSampleDist;
  String filename;
  List<_IrradiancePoint> irradiancePoints = [];
  BBox octreeBounds;
  _SubsurfaceOctreeNode octree;

  List<LightSampleOffsets> lightSampleOffsets;
  List<BSDFSampleOffsets> bsdfSampleOffsets;


  static DipoleSubsurfaceIntegrator Create(ParamSet params) {
    int maxDepth = params.findOneInt('maxdepth', 5);
    double maxError = params.findOneFloat('maxerror', 0.05);
    double minDist = params.findOneFloat('minsampledistance', 0.25);
    String pointsfile = params.findOneFilename('pointsfile', '');
    if (RenderOverrides.QuickRender()) {
      maxError *= 4.0;
      minDist *= 4.0;
    }
    return new DipoleSubsurfaceIntegrator(maxDepth, maxError, minDist,
                                          pointsfile);
  }
}

class _SubsurfaceOctreeNode {
  _SubsurfaceOctreeNode() {
    isLeaf = true;
    sumArea = 0.0;
    ips = children;
    for (int i = 0; i < 8; ++i) {
      ips[i] = null;
    }
  }

  void insert(BBox nodeBound, _IrradiancePoint ip) {
    Point pMid = (nodeBound.pMin + nodeBound.pMax) * 0.5;
    if (isLeaf) {
      // Add _IrradiancePoint_ to leaf octree node
      for (int i = 0; i < 8; ++i) {
        if (ips[i] == null) {
          ips[i] = ip;
          return;
        }
      }

      // Convert leaf node to interior node, redistribute points
      isLeaf = false;
      List<_IrradiancePoint> localIps = new List<_IrradiancePoint>(8);
      for (int i = 0; i < 8; ++i) {
        localIps[i] = ips[i];
        children[i] = null;
      }

      for (int i = 0; i < 8; ++i) {
        _IrradiancePoint ip = localIps[i];
        // Add _IrradiancePoint_ _ip_ to interior octree node
        int child = (ip.p.x > pMid.x ? 4 : 0) +
                    (ip.p.y > pMid.y ? 2 : 0) +
                    (ip.p.z > pMid.z ? 1 : 0);
        if (children[child] == null) {
          children[child] = new _SubsurfaceOctreeNode();
        }

        BBox childBound = Octree.OctreeChildBound(child, nodeBound, pMid);

        children[child].insert(childBound, ip);
      }
      // fall through to interior case to insert the new point...
    }

    // Add _IrradiancePoint_ _ip_ to interior octree node
    int child = (ip.p.x > pMid.x ? 4 : 0) +
                (ip.p.y > pMid.y ? 2 : 0) +
                (ip.p.z > pMid.z ? 1 : 0);

    if (children[child] == null) {
      children[child] = new _SubsurfaceOctreeNode();
    }

    BBox childBound = Octree.OctreeChildBound(child, nodeBound, pMid);

    children[child].insert(childBound, ip);
  }

  void initHierarchy() {
    if (isLeaf) {
      // Init _SubsurfaceOctreeNode_ leaf from _IrradiancePoint_s
      double sumWt = 0.0;
      int i;
      for (i = 0; i < 8; ++i) {
        if (ips[i] == null) {
          break;
        }
        double wt = ips[i].E.luminance();
        E += ips[i].E;
        p += ips[i].p * wt;
        sumWt += wt;
        sumArea += ips[i].area;
      }

      if (sumWt > 0.0) {
        p /= sumWt;
      }

      if (i != 0) {
        E /= i;
      }
    } else {
      // Init interior _SubsurfaceOctreeNode_
      double sumWt = 0.0;
      int nChildren = 0;
      for (int i = 0; i < 8; ++i) {
        if (children[i] == null) {
          continue;
        }
        ++nChildren;
        children[i].initHierarchy();
        double wt = children[i].E.luminance();
        E += children[i].E;
        p += children[i].p * wt;
        sumWt += wt;
        sumArea += children[i].sumArea;
      }

      if (sumWt > 0.0) {
        p /= sumWt;
      }

      E /= nChildren;
    }
  }

  Spectrum Mo(BBox nodeBound, Point pt, _DiffusionReflectance Rd,
              double maxError) {
    // Compute Mo at node if error is low enough
    double dw = sumArea / Vector.DistanceSquared(pt, p);
    if (dw < maxError && !nodeBound.inside(pt)) {
      Stats.SUBSURFACE_ADDED_INTERIOR_CONTRIBUTION(this);
      return Rd(Vector.DistanceSquared(pt, p)) * E * sumArea;
    }

    // Otherwise compute Mo from points in leaf or recursively visit children
    Spectrum Mo = new Spectrum(0.0);
    if (isLeaf) {
      // Accumulate Mo from leaf node
      for (int i = 0; i < 8; ++i) {
        if (ips[i] == null) {
          break;
        }
        Stats.SUBSURFACE_ADDED_POINT_CONTRIBUTION(ips[i]);
        Mo += Rd(Vector.DistanceSquared(pt, ips[i].p)) * ips[i].E * ips[i].area;
      }
    } else {
      // Recursively visit children nodes to compute Mo
      Point pMid = (nodeBound.pMin + nodeBound.pMax) * 0.5;
      for (int child = 0; child < 8; ++child) {
        if (children[child] == null) {
          continue;
        }
        BBox childBound = Octree.OctreeChildBound(child, nodeBound, pMid);
        Mo += children[child].Mo(childBound, pt, Rd, maxError);
      }
    }

    return Mo;
  }

  Point p = new Point(0.0);
  bool isLeaf;
  Spectrum E = new Spectrum(0.0);
  double sumArea;
  List children = new List(8);
  List ips;
}

class _DiffusionReflectance {
  _DiffusionReflectance(Spectrum sigma_a, Spectrum sigmap_s, double eta) {
    A = (1.0 + Fdr(eta)) / (1.0 - Fdr(eta));
    sigmap_t = sigma_a + sigmap_s;
    sigma_tr = (sigma_a * sigmap_t * 3.0).sqrt();
    alphap = sigmap_s / sigmap_t;
    zpos = Spectrum.ONE / sigmap_t;
    zneg = -zpos * (1.0 + (4.0 / 3.0) * A);
  }

  Spectrum call(double d2) {
    Spectrum dpos = (new Spectrum(d2) + zpos * zpos).sqrt();
    Spectrum dneg = (new Spectrum(d2) + zneg * zneg).sqrt();
    Spectrum Rd = (alphap / (4.0 * Math.PI)) *
                  ((zpos * (dpos * sigma_tr + Spectrum.ONE) *
                   (-sigma_tr * dpos).exp()) / (dpos * dpos * dpos) -
                   (zneg * (dneg * sigma_tr + Spectrum.ONE) *
                   (-sigma_tr * dneg).exp()) / (dneg * dneg * dneg));
    return Rd.clamp();
  }

  Spectrum zpos;
  Spectrum zneg;
  Spectrum sigmap_t;
  Spectrum sigma_tr;
  Spectrum alphap;
  double A;
}

class _IrradiancePoint {
  _IrradiancePoint(SurfacePoint sp, this.E) :
    p = sp.p,
    n = sp.n,
    area = sp.area,
    rayEpsilon = sp.rayEpsilon;

  Point p;
  Normal n;
  Spectrum E;
  double area;
  double rayEpsilon;
}
