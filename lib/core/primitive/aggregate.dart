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
 * An aggregate is a container for sub-primitives, and is used as the base
 * class for ray-trace acceleration structures.
 */
abstract class Aggregate extends Primitive {
  AreaLight getAreaLight() {
    LogSevere('Aggregate.getAreaLight() method'
             'called; should have gone to GeometricPrimitive');
    return null;
  }

  BSDF getBSDF(DifferentialGeometry dg, Transform xform)  {
    LogSevere('Aggregate.getBSDF() method'
              'called; should have gone to GeometricPrimitive');
    return null;
  }

  BSSRDF getBSSRDF(DifferentialGeometry dg, Transform xform) {
    LogSevere('Aggregate.getBSSRDF() method'
              'called; should have gone to GeometricPrimitive');
    return null;
  }
}
