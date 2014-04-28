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
 * Base class for cameras that use a projective transformation.
 */
abstract class ProjectiveCamera extends Camera {
  Transform cameraToScreen;
  Transform rasterToCamera;
  Transform screenToRaster;
  Transform rasterToScreen;
  double lensRadius;
  double focalDistance;

  ProjectiveCamera(AnimatedTransform cam2world,
                   this.cameraToScreen, List<double> screenWindow,
                   double sopen, double sclose, this.lensRadius,
                   this.focalDistance, Film film)
      : super(cam2world, sopen, sclose, film) {
    // Compute projective camera screen transformations
    screenToRaster =
         Transform.Scale(film.xResolution.toDouble(),
                         film.yResolution.toDouble(),
                         1.0) *
         Transform.Scale(1.0 / (screenWindow[1] - screenWindow[0]),
                         1.0 / (screenWindow[2] - screenWindow[3]),
                         1.0) *
         Transform.Translate(new Vector(-screenWindow[0].toDouble(),
                                        -screenWindow[3].toDouble(),
                                        0.0));

    rasterToScreen = Transform.Inverse(screenToRaster);
    rasterToCamera = Transform.Inverse(cameraToScreen) * rasterToScreen;
  }
}
