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
 * A sample of the camera image plane. This stores the sample's position on the
 * image plane, it's position on the camera lens, and the time as a value
 * between the opened and closed shutter times of the camera.
 */
class CameraSample {
  double imageX;
  double imageY;
  double lensU;
  double lensV;
  double time;

  CameraSample()
      : imageX = 0.0,
        imageY = 0.0,
        lensU = 0.0,
        lensV = 0.0,
        time = 0.0;

  CameraSample.from(CameraSample s)
      : imageX = s.imageX,
        imageY = s.imageY,
        lensU = s.lensU,
        lensV = s.lensV,
        time = s.time;
}
