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
part of cameras;

class OrthographicCamera extends ProjectiveCamera {
  OrthographicCamera(AnimatedTransform cam2world, List<double> screenWindow,
                     double sopen, double sclose, double lensr, double focald,
                     Film film)
      : super(cam2world, Transform.Orthographic(0.0, 1.0), screenWindow,
              sopen, sclose, lensr, focald, film) {
    // Compute differential changes in origin for ortho camera rays
    dxCamera = rasterToCamera.transformVector(new Vector(1.0, 0.0, 0.0));
    dyCamera = rasterToCamera.transformVector(new Vector(0.0, 1.0, 0.0));
  }

  double generateRay(CameraSample sample, Ray ray) {
    // Generate raster and camera samples
    Point Pras = new Point(sample.imageX, sample.imageY, 0.0);
    Point Pcamera = new Point();
    rasterToCamera.transformPoint(Pras, Pcamera);

    ray.set(Pcamera, new Vector(0.0, 0.0, 1.0), 0.0, INFINITY, sample.time);

    // Modify ray for depth of field
    if (lensRadius > 0.0) {
      // Sample point on lens
      List<double> lensU = [0.0];
      List<double> lensV = [0.0];
      ConcentricSampleDisk(sample.lensU, sample.lensV, lensU, lensV);
      lensU[0] *= lensRadius;
      lensV[0] *= lensRadius;

      // Compute point on plane of focus
      double ft = focalDistance / ray.direction.z;
      Point Pfocus = ray.pointAt(ft);

      // Update ray for effect of lens
      ray.origin = new Point(lensU[0], lensV[0], 0.0);
      ray.direction = Vector.Normalize(Pfocus - ray.origin);
    }

    cameraToWorld.transformRay(ray, ray);

    return 1.0;
  }

  double generateRayDifferential(CameraSample sample, RayDifferential ray) {
    // Compute main orthographic viewing ray

    // Generate raster and camera samples
    Point Pras = new Point(sample.imageX, sample.imageY, 0.0);
    Point Pcamera = new Point();
    rasterToCamera.transformPoint(Pras, Pcamera);
    ray.set(Pcamera, new Vector(0.0, 0.0, 1.0), 0.0, INFINITY, sample.time);

    // Modify ray for depth of field
    if (lensRadius > 0.0) {
      // Sample point on lens
      List<double> lensU = [0.0];
      List<double> lensV = [0.0];
      ConcentricSampleDisk(sample.lensU, sample.lensV, lensU, lensV);
      lensU[0] *= lensRadius;
      lensV[0] *= lensRadius;

      // Compute point on plane of focus
      double ft = focalDistance / ray.direction.z;
      Point Pfocus = ray.pointAt(ft);

      // Update ray for effect of lens
      ray.origin = new Point(lensU[0], lensV[0], 0.0);
      ray.direction = Vector.Normalize(Pfocus - ray.origin);
    }

    ray.rxOrigin = ray.origin + dxCamera;
    ray.ryOrigin = ray.origin + dyCamera;
    ray.rxDirection = ray.ryDirection = ray.direction;
    ray.hasDifferentials = true;
    cameraToWorld.transformRay(ray, ray);

    return 1.0;
  }

  static OrthographicCamera Create(ParamSet params, AnimatedTransform cam2world,
                                   Film film) {
    // Extract common camera parameters from _ParamSet_
    double shutteropen = params.findOneFloat("shutteropen", 0.0);
    double shutterclose = params.findOneFloat("shutterclose", 1.0);
    if (shutterclose < shutteropen) {
      LogWarning('Shutter close time [$shutterclose] < '
                 'shutter open [$shutteropen].  Swapping them.');
      double t = shutterclose;
      shutterclose = shutteropen;
      shutteropen = t;
    }

    double lensradius = params.findOneFloat("lensradius", 0.0);
    double focaldistance = params.findOneFloat("focaldistance", 1.0e30);
    double frame = params.findOneFloat("frameaspectratio",
                                       film.xResolution / film.yResolution);

    List<double> screen;

    List<double> sw = params.findFloat('screenwindow');
    if (sw != null && sw.length == 4) {
      screen = sw;
    } else {
      screen = [0.0, 0.0, 0.0, 0.0];
      if (frame > 1.0) {
        screen[0] = -frame;
        screen[1] =  frame;
        screen[2] = -1.0;
        screen[3] =  1.0;
      } else {
        screen[0] = -1.0;
        screen[1] =  1.0;
        screen[2] = -1.0 / frame;
        screen[3] =  1.0 / frame;
      }
    }

    return new OrthographicCamera(cam2world, screen, shutteropen, shutterclose,
                                  lensradius, focaldistance, film);
  }

  Vector dxCamera;
  Vector dyCamera;
}
