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
 * Base class for cameras, which generate rays for samples.
 */
abstract class Camera {
  AnimatedTransform cameraToWorld;
  double shutterOpen;
  double shutterClose;
  Film film;

  Camera(this.cameraToWorld, this.shutterOpen, this.shutterClose, this.film);

  double generateRay(CameraSample sample, Ray ray);

  double generateRayDifferential(CameraSample sample, RayDifferential rd) {
    double wt = generateRay(sample, rd);

    // Find ray after shifting one pixel in the $x$ direction
    CameraSample sshift = new CameraSample.from(sample);
    sshift.imageX++;
    Ray rx = new Ray();
    double wtx = generateRay(sshift, rx);
    rd.rxOrigin = new Point.from(rx.origin);
    rd.rxDirection = new Vector.from(rx.direction);

    // Find ray after shifting one pixel in the $y$ direction
    sshift.imageX--;
    sshift.imageY++;
    Ray ry = new Ray();
    double wty = generateRay(sshift, ry);
    rd.ryOrigin = new Point.from(ry.origin);
    rd.ryDirection = new Vector.from(ry.direction);
    if (wtx == 0.0 || wty == 0.0) {
      return 0.0;
    }

    rd.hasDifferentials = true;

    return wt;
  }
}
