/****************************************************************************
 *  Copyright (C) 2014 by Brendan Duncan.                                   *
 *                                                                          *
 *  This file is part of DartRay.                                           *
 *                                                                          *
 *  Licensed under the Apache License, Version 2.0 (the "License");         *
 *  you may not use this file except in compliance with the License.        *
 *  You may obtain a copy of the License at                                 *
 *                                                                          *
 *  http://www.apache.org/licenses/LICENSE-2.0                              *
 *                                                                          *
 *  Unless required by applicable law or agreed to in writing, software     *
 *  distributed under the License is distributed on an "AS IS" BASIS,       *
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 *  See the License for the specific language governing permissions and     *
 *  limitations under the License.                                          *
 *                                                                          *
 *   This project is based on PBRT v2 ; see http://www.pbrt.org             *
 *   pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.*
 ****************************************************************************/
part of core;

/**
 * Stores information about a point of intersection.
 */
class Intersection {
  Intersection() :
    dg = new DifferentialGeometry() {
    primitive = null;
    shapeId = 0;
    primitiveId = 0;
    rayEpsilon = 0.0;
  }

  void copy(Intersection other) {
    dg = other.dg;
    primitive = other.primitive;
    worldToObject = other.worldToObject;
    objectToWorld = other.objectToWorld;
    shapeId = other.shapeId;
    primitiveId = other.primitiveId;
    rayEpsilon = other.rayEpsilon;
  }

  BSDF getBSDF(RayDifferential ray) {
    dg.computeDifferentials(ray);
    BSDF bsdf = primitive.getBSDF(dg, objectToWorld);
    return bsdf;
  }

  BSSRDF getBSSRDF(RayDifferential ray) {
    dg.computeDifferentials(ray);
    BSSRDF bssrdf = primitive.getBSSRDF(dg, objectToWorld);
    return bssrdf;
  }

  Spectrum Le(Vector wo) {
    AreaLight area = primitive.getAreaLight();
    return area != null ? area.L(dg.p, dg.nn, wo) : new Spectrum(0.0);
  }

  DifferentialGeometry dg;
  Primitive primitive;
  Transform worldToObject;
  Transform objectToWorld;
  int shapeId;
  int primitiveId;
  double rayEpsilon;
}
