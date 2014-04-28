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
 * Base class for Material extensions.
 */
abstract class Material {
  BSDF getBSDF(DifferentialGeometry dgGeom,
               DifferentialGeometry dgShading);

  BSSRDF getBSSRDF(DifferentialGeometry dgGeom,
                   DifferentialGeometry dgShading) {
    return null;
  }

  static void Bump(Texture d, DifferentialGeometry dgGeom,
                   DifferentialGeometry dgs,
                   DifferentialGeometry dgBump) {
    // Compute offset positions and evaluate displacement texture
    DifferentialGeometry dgEval = new DifferentialGeometry.from(dgs);

    // Shift dgEval du in the u direction
    double du = 0.5 * (dgs.dudx.abs() + dgs.dudy.abs());
    if (du == 0.0) {
      du = 0.01;
    }

    dgEval.p = dgs.p + dgs.dpdu * du;
    dgEval.u = dgs.u + du;
    dgEval.nn = Normal.Normalize(Vector.Cross(dgs.dpdu, dgs.dpdv) +
                                 dgs.dndu * du);
    double uDisplace = d.evaluate(dgEval);

    // Shift dgEval dv in the v direction
    double dv = 0.5 * (dgs.dvdx.abs() + dgs.dvdy.abs());
    if (dv == 0.0) {
      dv = 0.01;
    }
    dgEval.p = dgs.p + dgs.dpdv * dv;
    dgEval.u = dgs.u;
    dgEval.v = dgs.v + dv;
    dgEval.nn = Normal.Normalize(Vector.Cross(dgs.dpdu, dgs.dpdv) +
                                 dgs.dndv * dv);
    double vDisplace = d.evaluate(dgEval);
    double displace = d.evaluate(dgs);

    // Compute bump-mapped differential geometry
    dgBump.copy(dgs);

    dgBump.dpdu = dgs.dpdu +
                  new Vector.from(dgs.nn) * (uDisplace - displace) / du +
                  new Vector.from(dgs.dndu) * displace;

    dgBump.dpdv = dgs.dpdv +
                  new Vector.from(dgs.nn) * (vDisplace - displace) / dv +
                  new Vector.from(dgs.dndv) * displace;

    dgBump.nn = new Normal.from(Vector.Normalize(Vector.Cross(dgBump.dpdu,
                                                              dgBump.dpdv)));

    // dgs.shape.reverseOrientation ^ dgs.shape.transformSwapsHandedness
    if ((dgs.shape.reverseOrientation || dgs.shape.transformSwapsHandedness) &&
        (dgs.shape.reverseOrientation != dgs.shape.transformSwapsHandedness)) {
       dgBump.nn *= -1.0;
    }

    // Orient shading normal to match geometric normal
    dgBump.nn = Normal.FaceForward(dgBump.nn, dgGeom.nn);
  }
}
