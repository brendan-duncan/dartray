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

class PhotonMapIntegrator extends SurfaceIntegrator {
  PhotonMapIntegrator(this.nCausticPhotonsWanted, this.nIndirectPhotonsWanted,
                      this.nLookup, this.maxSpecularDepth,
                      this.maxPhotonDepth, double maxdist,
                      this.finalGather, this.gatherSamples,
                      double ga) :
    maxDistSquared = [maxdist * maxdist],
    cosGatherAngle = Math.cos(Radians(ga));

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

    L += Integrator.UniformSampleAllLights(scene, renderer, p, n,
                                           wo, isect.rayEpsilon, ray.time,
                                           bsdf, sample, rng,
                                           lightSampleOffsets,
                                           bsdfSampleOffsets);

    // Compute caustic lighting for photon map integrator
    List<ClosePhoton> lookupBuf = new List<ClosePhoton>(nLookup);

    L += _LPhoton(causticMap, nCausticPaths, nLookup, lookupBuf, bsdf,
                  rng, isect, wo, maxDistSquared);

    // Compute indirect lighting for photon map integrator
    if (finalGather && indirectMap != null) {
      // Do one-bounce final gather for photon map
      int nonSpecular = BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_DIFFUSE |
                        BSDF_GLOSSY;
      if (bsdf.numComponents(nonSpecular) > 0) {
        // Find indirect photons around point for importance sampling
        const int nIndirSamplePhotons = 50;
        PhotonProcess proc = new PhotonProcess(nIndirSamplePhotons,
                                    new List<ClosePhoton>(nIndirSamplePhotons));
        double searchDist2 = maxDistSquared[0];
        while (proc.nFound < nIndirSamplePhotons) {
          List<double> md2 = [searchDist2];
          proc.nFound = 0;
          indirectMap.lookup(p, proc, md2);
          searchDist2 *= 2.0;
        }

        // Copy photon directions to local array
        List<Vector> photonDirs = new List<Vector>(nIndirSamplePhotons);
        for (int i = 0; i < nIndirSamplePhotons; ++i) {
          photonDirs[i] = proc.photons[i].photon.wi;
        }

        // Use BSDF to do final gathering
        Spectrum Li = new Spectrum(0.0);
        for (int i = 0; i < gatherSamples; ++i) {
          // Sample random direction from BSDF for final gather ray
          Vector wi = new Vector();
          List<double> pdf = [0.0];
          BSDFSample bsdfSample = new BSDFSample.sample(sample,
                                                    bsdfGatherSampleOffsets, i);
          Spectrum fr = bsdf.sample_f(wo, wi, bsdfSample,
                                      pdf, BSDF_ALL & ~BSDF_SPECULAR);
          if (fr.isBlack() || pdf[0] == 0.0) {
            continue;
          }
          assert(pdf[0] >= 0.0);

          // Trace BSDF final gather ray and accumulate radiance
          RayDifferential bounceRay = new RayDifferential.child(p, wi, ray,
                                                              isect.rayEpsilon);
          Intersection gatherIsect = new Intersection();
          if (scene.intersect(bounceRay, gatherIsect)) {
            // Compute exitant radiance _Lindir_ using radiance photons
            Spectrum Lindir = new Spectrum(0.0);
            Normal nGather = gatherIsect.dg.nn;
            nGather = Normal.FaceForward(nGather, -bounceRay.direction);
            RadiancePhotonProcess proc = new RadiancePhotonProcess(nGather);
            List<double> md2 = [INFINITY];
            radianceMap.lookup(gatherIsect.dg.p, proc, md2);
            if (proc.photon != null) {
              Lindir = proc.photon.Lo;
            }
            Lindir *= renderer.transmittance(scene, bounceRay, null, rng);

            // Compute MIS weight for BSDF-sampled gather ray

            // Compute PDF for photon-sampling of direction _wi_
            double photonPdf = 0.0;
            double conePdf = UniformConePdf(cosGatherAngle);
            for (int j = 0; j < nIndirSamplePhotons; ++j) {
              if (Vector.Dot(photonDirs[j], wi) > 0.999 * cosGatherAngle) {
                photonPdf += conePdf;
              }
            }
            photonPdf /= nIndirSamplePhotons;
            double wt = PowerHeuristic(gatherSamples, pdf[0], gatherSamples,
                                       photonPdf);

            Li += fr * Lindir * (Vector.AbsDot(wi, n) * wt / pdf[0]);
          }
        }

        L += Li / gatherSamples;

        // Use nearby photons to do final gathering
        Li = new Spectrum(0.0);
        for (int i = 0; i < gatherSamples; ++i) {
          // Sample random direction using photons for final gather ray
          BSDFSample gatherSample =
                     new BSDFSample.sample(sample, indirGatherSampleOffsets, i);
          int photonNum = Math.min(nIndirSamplePhotons - 1,
                     (gatherSample.uComponent * nIndirSamplePhotons).floor());

          // Sample gather ray direction from _photonNum_
          Vector vx = new Vector();
          Vector vy = new Vector();
          Vector.CoordinateSystem(photonDirs[photonNum], vx, vy);
          Vector wi = UniformSampleCone2(gatherSample.uDir[0], gatherSample.uDir[1],
                                        cosGatherAngle, vx, vy,
                                        photonDirs[photonNum]);

          // Trace photon-sampled final gather ray and accumulate radiance
          Spectrum fr = bsdf.f(wo, wi);
          if (fr.isBlack()) {
            continue;
          }

          RayDifferential bounceRay = new RayDifferential.child(p, wi, ray,
                                                              isect.rayEpsilon);
          Intersection gatherIsect = new Intersection();
          Stats.PHOTON_MAP_STARTED_GATHER_RAY(bounceRay);
          if (scene.intersect(bounceRay, gatherIsect)) {
            // Compute exitant radiance _Lindir_ using radiance photons
            Spectrum Lindir = new Spectrum(0.0);
            Normal nGather = gatherIsect.dg.nn;
            nGather = Normal.FaceForward(nGather, -bounceRay.direction);
            RadiancePhotonProcess proc = new RadiancePhotonProcess(nGather);
            List<double> md2 = [INFINITY];
            radianceMap.lookup(gatherIsect.dg.p, proc, md2);
            if (proc.photon != null) {
              Lindir = proc.photon.Lo;
            }
            Lindir *= renderer.transmittance(scene, bounceRay, null, rng);

            // Compute PDF for photon-sampling of direction _wi_
            double photonPdf = 0.0;
            double conePdf = UniformConePdf(cosGatherAngle);
            for (int j = 0; j < nIndirSamplePhotons; ++j) {
              if (Vector.Dot(photonDirs[j], wi) > 0.999 * cosGatherAngle) {
                photonPdf += conePdf;
              }
            }
            photonPdf /= nIndirSamplePhotons;

            // Compute MIS weight for photon-sampled gather ray
            double bsdfPdf = bsdf.pdf(wo, wi);
            double wt = PowerHeuristic(gatherSamples, photonPdf, gatherSamples,
                                       bsdfPdf);
            Li += fr * Lindir * Vector.AbsDot(wi, n) * wt / photonPdf;
          }
          Stats.PHOTON_MAP_FINISHED_GATHER_RAY(bounceRay);
        }

        L += Li / gatherSamples;
      }
    } else {
      L += _LPhoton(indirectMap, nIndirectPaths, nLookup, lookupBuf,
                    bsdf, rng, isect, wo, maxDistSquared);
    }

    if (ray.depth + 1 < maxSpecularDepth) {
      Vector wi = new Vector();
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

    // Request samples for final gathering
    if (finalGather) {
      gatherSamples = Math.max(1, gatherSamples ~/ 2);
      if (sampler != null) {
        gatherSamples = sampler.roundSize(gatherSamples);
      }
      bsdfGatherSampleOffsets = new BSDFSampleOffsets(gatherSamples, sample);
      indirGatherSampleOffsets = new BSDFSampleOffsets(gatherSamples, sample);
    }
  }

  void preprocess(Scene scene, Camera camera, Renderer renderer) {
    if (scene.lights.length == 0) {
      return;
    }

    // Declare shared variables for photon shooting
    List<int> nDirectPaths = [0];
    List<Photon> causticPhotons = [];
    List<Photon> directPhotons = [];
    List<Photon> indirectPhotons = [];
    List<RadiancePhoton> radiancePhotons = [];
    List<bool> abortTasks = [false];

    List<int> nshot = [0];
    List<Spectrum> rpReflectances = [];
    List<Spectrum> rpTransmittances = [];

    // Compute light power CDF for photon shooting
    Distribution1D lightDistribution = Integrator.ComputeLightSamplingCDF(scene);

    PhotonShootingTask task = new PhotonShootingTask(
                0, camera != null ? camera.shutterOpen : 0.0,
                this, abortTasks, nDirectPaths,
                directPhotons, indirectPhotons, causticPhotons, radiancePhotons,
                rpReflectances, rpTransmittances,
                nshot, lightDistribution, scene, renderer);
    task.run();

    // Build kd-trees for indirect and caustic photons
    KdTree directMap;
    if (directPhotons.isNotEmpty) {
      directMap = new KdTree(directPhotons);
    }
    if (causticPhotons.isNotEmpty) {
      causticMap = new KdTree(causticPhotons);
    }
    if (indirectPhotons.isNotEmpty) {
      indirectMap = new KdTree(indirectPhotons);
    }

    // Precompute radiance at a subset of the photons
    if (finalGather && radiancePhotons.isNotEmpty) {
      // Launch tasks to compute photon radiances
      ComputeRadianceTask task = new ComputeRadianceTask(0, 1, radiancePhotons,
              rpReflectances, rpTransmittances,
              nLookup, maxDistSquared, nDirectPaths[0], directMap,
              nIndirectPaths, indirectMap,
              nCausticPaths, causticMap);
      task.run();

      radianceMap = new KdTree(radiancePhotons);
    }
  }

  static PhotonMapIntegrator Create(ParamSet params) {
    int nCaustic = params.findOneInt('causticphotons', 20000);
    int nIndirect = params.findOneInt('indirectphotons', 100000);
    int nUsed = params.findOneInt('nused', 50);

    int maxSpecularDepth = params.findOneInt('maxspeculardepth', 5);
    int maxPhotonDepth = params.findOneInt('maxphotondepth', 5);
    bool finalGather = params.findOneBool('finalgather', true);
    int gatherSamples = params.findOneInt('finalgathersamples', 32);

    double maxDist = params.findOneFloat('maxdist', 0.1);
    double gatherAngle = params.findOneFloat('gatherangle', 10.0);

    return new PhotonMapIntegrator(nCaustic, nIndirect, nUsed,
                                   maxSpecularDepth, maxPhotonDepth, maxDist,
                                   finalGather, gatherSamples, gatherAngle);
  }

  int nCausticPhotonsWanted;
  int nIndirectPhotonsWanted;
  int nLookup;
  List<double> maxDistSquared;
  int maxSpecularDepth;
  int maxPhotonDepth;
  bool finalGather;
  int gatherSamples;
  double cosGatherAngle;

  // Declare sample parameters for light source sampling
  List<LightSampleOffsets> lightSampleOffsets;
  List<BSDFSampleOffsets> bsdfSampleOffsets;
  BSDFSampleOffsets bsdfGatherSampleOffsets;
  BSDFSampleOffsets indirGatherSampleOffsets;
  int nCausticPaths = 0;
  int nIndirectPaths = 0;
  KdTree causticMap;
  KdTree indirectMap;
  KdTree radianceMap;
}

class Photon {
  Photon(this.p, this.alpha, this.wi);

  Point p;
  Spectrum alpha;
  Vector wi;
}

class RadiancePhoton {
  RadiancePhoton(this.p, this.n);

  Point p;
  Normal n;
  Spectrum Lo = new Spectrum(0.0);
}

class PhotonShootingTask {
  PhotonShootingTask(this.taskNum, this.time, this.integrator,
      this.abortTasks, this.nDirectPaths,
      this.directPhotons, this.indirectPhotons, this.causticPhotons,
      this.radiancePhotons, this.rpReflectances, this.rpTransmittances,
      this.nshot, this.lightDistribution, this.scene,
      this.renderer);

  void run() {
    // Declare local variables for _PhotonShootingTask_
    RNG rng = new RNG(31 * taskNum);
    List<Photon> localDirectPhotons = [];
    List<Photon> localIndirectPhotons = [];
    List<Photon> localCausticPhotons = [];
    List<RadiancePhoton> localRadiancePhotons = [];
    int totalPaths = 0;
    bool causticDone = (integrator.nCausticPhotonsWanted == 0);
    bool indirectDone = (integrator.nIndirectPhotonsWanted == 0);

    PermutedHalton halton = new PermutedHalton(6, rng);
    List<Spectrum> localRpReflectances = [];
    List<Spectrum> localRpTransmittances = [];
    List<double> u = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    List<double> lightPdf = [0.0];

    while (true) {
      // Follow photon paths for a block of samples
      const int blockSize = 4096;
      for (int i = 0; i < blockSize; ++i) {
        halton.sample(++totalPaths, u);
        // Choose light to shoot photon from
        int lightNum = lightDistribution.sampleDiscrete(u[0], lightPdf);
        Light light = scene.lights[lightNum];

        // Generate _photonRay_ from light source and initialize _alpha_
        RayDifferential photonRay = new RayDifferential();
        List<double> pdf = [0.0];
        LightSample ls = new LightSample(u[1], u[2], u[3]);
        Normal Nl = new Normal();

        Spectrum Le = light.sampleL(scene, ls, u[4], u[5], time, photonRay,
                                     Nl, pdf);

        if (pdf[0] == 0.0 || Le.isBlack()) {
          continue;
        }

        Spectrum alpha = (Le * Vector.AbsDot(Nl, photonRay.direction)) /
                         (pdf[0] * lightPdf[0]);
        if (!alpha.isBlack()) {
          // Follow photon path through scene and record intersections
          Stats.PHOTON_MAP_STARTED_RAY_PATH(photonRay, alpha);
          bool specularPath = true;
          Intersection photonIsect = new Intersection();
          int nIntersections = 0;
          while (scene.intersect(photonRay, photonIsect)) {
            ++nIntersections;
            // Handle photon/surface intersection
            alpha *= renderer.transmittance(scene, photonRay, null, rng);
            BSDF photonBSDF = photonIsect.getBSDF(photonRay);
            int specularType = BSDF_REFLECTION | BSDF_TRANSMISSION |
                               BSDF_SPECULAR;
            bool hasNonSpecular = (photonBSDF.numComponents() >
                                   photonBSDF.numComponents(specularType));
            Vector wo = -photonRay.direction;
            if (hasNonSpecular) {
              // Deposit photon at surface
              Photon photon = new Photon(photonIsect.dg.p, alpha, wo);
              bool depositedPhoton = false;
              if (specularPath && nIntersections > 1) {
                if (!causticDone) {
                  Stats.PHOTON_MAP_DEPOSITED_CAUSTIC_PHOTON(photonIsect.dg,
                                                            alpha, wo);
                  depositedPhoton = true;
                  localCausticPhotons.add(photon);
                }
              } else {
                // Deposit either direct or indirect photon
                // stop depositing direct photons once indirectDone is true;
                // don't want to waste memory storing too many if we're going a
                // long time trying to get enough caustic photons desposited.
                if (nIntersections == 1 && !indirectDone &&
                    integrator.finalGather) {
                  Stats.PHOTON_MAP_DEPOSITED_DIRECT_PHOTON(photonIsect.dg,
                                                           alpha, wo);
                  depositedPhoton = true;
                  localDirectPhotons.add(photon);
                } else if (nIntersections > 1 && !indirectDone) {
                  Stats.PHOTON_MAP_DEPOSITED_INDIRECT_PHOTON(photonIsect.dg,
                                                             alpha, wo);
                  depositedPhoton = true;
                  localIndirectPhotons.add(photon);
                }
              }

              // Possibly create radiance photon at photon intersection point
              if (depositedPhoton && integrator.finalGather &&
                  rng.randomFloat() < 0.125) {
                Normal n = photonIsect.dg.nn;
                n = Normal.FaceForward(n, -photonRay.direction);
                localRadiancePhotons.add(new RadiancePhoton(photonIsect.dg.p, n));
                Spectrum rho_r = photonBSDF.rho(rng, BSDF_ALL_REFLECTION);
                localRpReflectances.add(rho_r);
                Spectrum rho_t = photonBSDF.rho(rng, BSDF_ALL_TRANSMISSION);
                localRpTransmittances.add(rho_t);
              }
            }

            if (nIntersections >= integrator.maxPhotonDepth) {
              break;
            }

            // Sample new photon ray direction
            Vector wi = new Vector();
            List<double> pdf = [0.0];
            List<int> flags = [0];
            Spectrum fr = photonBSDF.sample_f(wo, wi, new BSDFSample.random(rng),
                                              pdf, BSDF_ALL, flags);

            if (fr.isBlack() || pdf[0] == 0.0) {
              break;
            }

            Spectrum anew = alpha * fr *
                            Vector.AbsDot(wi, photonBSDF.dgShading.nn) / pdf[0];

            // Possibly terminate photon path with Russian roulette
            double continueProb = Math.min(1.0,
                                          anew.luminance() / alpha.luminance());
            if (rng.randomFloat() > continueProb) {
              break;
            }

            alpha = anew / continueProb;
            specularPath = specularPath && ((flags[0] & BSDF_SPECULAR) != 0);

            if (indirectDone && !specularPath) {
              break;
            }

            photonRay = new RayDifferential.child(photonIsect.dg.p, wi,
                                                       photonRay,
                                                       photonIsect.rayEpsilon);
          }
          Stats.PHOTON_MAP_FINISHED_RAY_PATH(photonRay, alpha);
        }
      }

      // Merge local photon data with data in _PhotonIntegrator_
      {
        // Give up if we're not storing enough photons
        if (abortTasks[0]) {
          return;
        }

        if (nshot[0] > 500000 &&
            (_unsuccessful(integrator.nCausticPhotonsWanted,
                           causticPhotons.length, blockSize) ||
             _unsuccessful(integrator.nIndirectPhotonsWanted,
                           indirectPhotons.length, blockSize))) {
          LogError('Unable to store enough photons. Giving up.\n');
          causticPhotons.clear();
          indirectPhotons.clear();
          radiancePhotons.clear();
          abortTasks[0] = true;
          return;
        }

        nshot[0] += blockSize;

        // Merge indirect photons into shared array
        if (!indirectDone) {
          integrator.nIndirectPaths += blockSize;
          for (int i = 0; i < localIndirectPhotons.length; ++i) {
            indirectPhotons.add(localIndirectPhotons[i]);
          }
          localIndirectPhotons.clear();

          if (indirectPhotons.length >= integrator.nIndirectPhotonsWanted) {
            indirectDone = true;
          }

          nDirectPaths[0] += blockSize;
          for (int i = 0; i < localDirectPhotons.length; ++i) {
            directPhotons.add(localDirectPhotons[i]);
          }
          localDirectPhotons.clear();
        }

        // Merge direct, caustic, and radiance photons into shared array
        if (!causticDone) {
          integrator.nCausticPaths += blockSize;
          for (int i = 0; i < localCausticPhotons.length; ++i) {
            causticPhotons.add(localCausticPhotons[i]);
          }
          localCausticPhotons.clear();
          if (causticPhotons.length >= integrator.nCausticPhotonsWanted) {
            causticDone = true;
          }
        }

        for (int i = 0; i < localRadiancePhotons.length; ++i) {
          radiancePhotons.add(localRadiancePhotons[i]);
        }
        localRadiancePhotons.clear();

        for (int i = 0; i < localRpReflectances.length; ++i) {
          rpReflectances.add(localRpReflectances[i]);
        }
        localRpReflectances.clear();
        for (int i = 0; i < localRpTransmittances.length; ++i) {
          rpTransmittances.add(localRpTransmittances[i]);
        }
        localRpTransmittances.clear();
      }

      // Exit task if enough photons have been found
      if (indirectDone && causticDone) {
        break;
      }
    }
  }

  int taskNum;
  double time;
  PhotonMapIntegrator integrator;
  List<bool> abortTasks;
  List<int> nDirectPaths;
  List<Photon> directPhotons;
  List<Photon> indirectPhotons;
  List<Photon> causticPhotons;
  List<RadiancePhoton> radiancePhotons;
  List<Spectrum> rpReflectances;
  List<Spectrum> rpTransmittances;
  List<int> nshot;
  Distribution1D lightDistribution;
  Scene scene;
  Renderer renderer;
}

class ComputeRadianceTask {
  ComputeRadianceTask(this.taskNum, this.numTasks,
  this.radiancePhotons, this.rpReflectances,
  this.rpTransmittances,
  this.nLookup, this.maxDistSquared,
  this.nDirectPaths, this.directMap,
  this.nIndirectPaths, this.indirectMap,
  this.nCausticPaths, this.causticMap);

  void run() {
    // Compute range of radiance photons to process in task
    int taskSize = radiancePhotons.length ~/ numTasks;
    int excess = radiancePhotons.length % numTasks;
    int rpStart = Math.min(taskNum, excess) * (taskSize + 1) +
                  Math.max(0, taskNum - excess) * taskSize;
    int rpEnd = rpStart + taskSize + (taskNum < excess ? 1 : 0);

    if (taskNum == numTasks - 1) {
      assert(rpEnd == radiancePhotons.length);
    }

    List<ClosePhoton> lookupBuf = new List<ClosePhoton>(nLookup);
    for (int i = rpStart; i < rpEnd; ++i) {
      // Compute radiance for radiance photon _i_
      RadiancePhoton rp = radiancePhotons[i];
      Spectrum rho_r = rpReflectances[i];
      Spectrum rho_t = rpTransmittances[i];
      if (!rho_r.isBlack()) {
        // Accumulate outgoing radiance due to reflected irradiance
        Spectrum E = _EPhoton(directMap, nDirectPaths, nLookup, lookupBuf,
                              maxDistSquared, rp.p, rp.n) +
                     _EPhoton(indirectMap, nIndirectPaths, nLookup, lookupBuf,
                              maxDistSquared, rp.p, rp.n) +
                     _EPhoton(causticMap, nCausticPaths, nLookup, lookupBuf,
                              maxDistSquared, rp.p, rp.n);
        rp.Lo += rho_r * E * INV_PI;
      }

      if (!rho_t.isBlack()) {
        // Accumulate outgoing radiance due to transmitted irradiance
        Spectrum E = _EPhoton(directMap, nDirectPaths, nLookup, lookupBuf,
                              maxDistSquared, rp.p, -rp.n) +
                     _EPhoton(indirectMap, nIndirectPaths, nLookup, lookupBuf,
                              maxDistSquared, rp.p, -rp.n) +
                     _EPhoton(causticMap, nCausticPaths, nLookup, lookupBuf,
                              maxDistSquared, rp.p, -rp.n);
        rp.Lo += rho_t * E * INV_PI;
      }
    }
  }

  int taskNum;
  int numTasks;
  List<RadiancePhoton> radiancePhotons;
  List<Spectrum> rpReflectances;
  List<Spectrum> rpTransmittances;
  int nLookup;
  List<double> maxDistSquared;
  int nDirectPaths;
  int nIndirectPaths;
  int nCausticPaths;
  KdTree directMap;
  KdTree indirectMap;
  KdTree causticMap;
}

class ClosePhoton {
  ClosePhoton([this.photon, this.distanceSquared]);

  bool operator<(ClosePhoton p2) {
    return distanceSquared == p2.distanceSquared ?
        (photon.hashCode < p2.photon.hashCode) :
        (distanceSquared < p2.distanceSquared);
  }

  Photon photon;
  double distanceSquared;
}

class PhotonProcess {
  PhotonProcess(this.nLookup, this.photons) :
    nFound = 0;

  void call(Point p, Photon photon, double distSquared,
            List<double> maxDistSquared) {
    if (nFound < nLookup) {
      // Add photon to unordered array of photons
      photons[nFound++] = new ClosePhoton(photon, distSquared);
      if (nFound == nLookup) {
        make_heap(photons, 0, nLookup);
        maxDistSquared[0] = photons[0].distanceSquared;
      }
    } else {
      // Remove most distant photon from heap and add new photon
      pop_heap(photons, 0, nLookup);
      photons[nLookup - 1] = new ClosePhoton(photon, distSquared);
      push_heap(photons, 0, nLookup);
      maxDistSquared[0] = photons[0].distanceSquared;
    }
  }

  List<ClosePhoton> photons;
  int nLookup;
  int nFound;
}

class RadiancePhotonProcess {
  RadiancePhotonProcess(this.n);

  void call(Point p, RadiancePhoton rp,
            double distSquared, List<double> maxDistSquared) {
    if (Vector.Dot(rp.n, n) > 0) {
      photon = rp;
      maxDistSquared[0] = distSquared;
    }
  }

  Normal n;
  RadiancePhoton photon;
}

bool _unsuccessful(int needed, int found, int shot) {
  return (found < needed && (found == 0 || found < shot / 1024));
}

double _kernel(Photon photon, Point p, double maxDist2) {
  double s = (1.0 - Vector.DistanceSquared(photon.p, p) / maxDist2);
  return 3.0 * INV_PI * s * s;
}

Spectrum _LPhoton(KdTree map, int nPaths, int nLookup,
      List<ClosePhoton> lookupBuf, BSDF bsdf, RNG rng,
      Intersection isect, Vector wo, List<double> maxDist2) {
  Spectrum L = new Spectrum(0.0);
  int nonSpecular = BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_DIFFUSE |
                    BSDF_GLOSSY;

  if (map != null && bsdf.numComponents(nonSpecular) > 0) {
    Stats.PHOTON_MAP_STARTED_LOOKUP(isect.dg);
    // Do photon map lookup at intersection point
    PhotonProcess proc = new PhotonProcess(nLookup, lookupBuf);
    map.lookup(isect.dg.p, proc, maxDist2);

    // Estimate reflected radiance due to incident photons
    List<ClosePhoton> photons = proc.photons;
    int nFound = proc.nFound;
    Normal Nf = Normal.FaceForward(bsdf.dgShading.nn, wo);
    if (bsdf.numComponents(BSDF_REFLECTION | BSDF_TRANSMISSION |
                           BSDF_GLOSSY) > 0) {
      // Compute exitant radiance from photons for glossy surface
      for (int i = 0; i < nFound; ++i) {
        Photon p = photons[i].photon;
        double k = _kernel(p, isect.dg.p, maxDist2[0]);
        L += bsdf.f(wo, p.wi) * p.alpha * (k / (nPaths * maxDist2[0]));
      }
    } else {
      // Compute exitant radiance from photons for diffuse surface
      Spectrum Lr = new Spectrum(0.0);
      Spectrum Lt = new Spectrum(0.0);
      for (int i = 0; i < nFound; ++i) {
        if (Vector.Dot(Nf, photons[i].photon.wi) > 0.0) {
          double k = _kernel(photons[i].photon, isect.dg.p, maxDist2[0]);
          Lr += photons[i].photon.alpha * (k / (nPaths * maxDist2[0]));
        } else {
          double k = _kernel(photons[i].photon, isect.dg.p, maxDist2[0]);
          Lt += photons[i].photon.alpha * (k / (nPaths * maxDist2[0]));
        }
      }
      L += Lr * bsdf.rho2(wo, rng, BSDF_ALL_REFLECTION) * INV_PI +
           Lt * bsdf.rho2(wo, rng, BSDF_ALL_TRANSMISSION) * INV_PI;
    }
    Stats.PHOTON_MAP_FINISHED_LOOKUP(isect.dg, proc.nFound, proc.nLookup, L);
  }
  return L;
}


Spectrum _EPhoton(KdTree map, int count, int nLookup,
                  List<ClosePhoton> lookupBuf, List<double> maxDist2, Point p,
                  Normal n) {
  if (map == null) {
    return new Spectrum(0.0);
  }

  // Lookup nearby photons at irradiance computation point
  PhotonProcess proc = new PhotonProcess(nLookup, lookupBuf);
  List<double> md2 = [maxDist2[0]];
  map.lookup(p, proc, md2);
  assert(md2[0] > 0.0);

  // Accumulate irradiance value from nearby photons
  if (proc.nFound == 0) {
    return new Spectrum(0.0);
  }

  List<ClosePhoton> photons = proc.photons;
  Spectrum E = new Spectrum(0.0);
  for (int i = 0; i < proc.nFound; ++i) {
    if (Vector.Dot(n, photons[i].photon.wi) > 0.0) {
      E += photons[i].photon.alpha;
    }
  }

  return E / (count * md2[0] * Math.PI);
}
