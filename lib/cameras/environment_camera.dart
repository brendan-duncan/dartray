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
part of cameras;

class EnvironmentCamera extends Camera {
  EnvironmentCamera(AnimatedTransform cam2world, double sopen,
                    double sclose, Film film)
      : super(cam2world, sopen, sclose, film);

  double generateRay(CameraSample sample, Ray ray) {
    // Compute environment camera ray direction
    double theta = Math.PI * sample.imageY / film.yResolution;
    double phi = 2 * Math.PI * sample.imageX / film.xResolution;
    Vector dir = new Vector(Math.sin(theta) * Math.cos(phi),
                            Math.cos(theta),
                            Math.sin(theta) * Math.sin(phi));
    ray.set(new Point(), dir, 0.0, INFINITY, sample.time);
    cameraToWorld.transformRay(ray, ray);
    return 1.0;
  }

  static EnvironmentCamera Create(ParamSet params, AnimatedTransform cam2world,
                                  Film film) {
    // Extract common camera parameters from _ParamSet_
    double shutteropen = params.findOneFloat('shutteropen', 0.0);
    double shutterclose = params.findOneFloat('shutterclose', 1.0);

    if (shutterclose < shutteropen) {
        LogWarning('Shutter close time [$shutterclose] < shutter open '
                   '[$shutteropen]. Swapping them.');
        double t = shutteropen;
        shutteropen = shutterclose;
        shutterclose = t;
    }

    double frame = params.findOneFloat('frameaspectratio',
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

    return new EnvironmentCamera(cam2world, shutteropen, shutterclose, film);
  }
}
