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
 * A primitive with an animatable world-space transformation.
 */
class TransformedPrimitive extends Primitive {
  TransformedPrimitive(this.primitive, AnimatedTransform w2p)
      : worldToPrimitive = new AnimatedTransform.from(w2p);

  bool intersect(Ray r, Intersection isect) {
    Transform w2p = new Transform();
    worldToPrimitive.interpolate(r.time, w2p);

    Ray ray = w2p.transformRay(r);
    if (!primitive.intersect(ray, isect)) {
      return false;
    }

    r.maxDistance = ray.maxDistance;
    isect.primitiveId = primitiveId;

    if (!w2p.isIdentity()) {
      // Compute world-to-object transformation for instance
      isect.worldToObject = isect.worldToObject * w2p;
      isect.objectToWorld = Transform.Inverse(isect.worldToObject);

      // Transform instance's differential geometry to world space
      Transform p2w = Transform.Inverse(w2p);
      isect.dg.p = p2w.transformPoint(isect.dg.p);
      isect.dg.nn = Vector.Normalize(p2w.transformNormal(isect.dg.nn));
      isect.dg.dpdu = p2w.transformVector(isect.dg.dpdu);
      isect.dg.dpdv = p2w.transformVector(isect.dg.dpdv);
      isect.dg.dndu = p2w.transformNormal(isect.dg.dndu);
      isect.dg.dndv = p2w.transformNormal(isect.dg.dndv);
    }

    return true;
  }

  bool intersectP(Ray r) {
    return primitive.intersectP(worldToPrimitive.transformRay(r));
  }

  AreaLight getAreaLight() {
    return null;
  }

  BSDF getBSDF(DifferentialGeometry dg, Transform ObjectToWorld) {
    return null;
  }

  BSSRDF getBSSRDF(DifferentialGeometry dg, Transform ObjectToWorld) {
    return null;
  }

  BBox worldBound() {
    return worldToPrimitive.motionBounds(primitive.worldBound(), true);
  }

  Primitive primitive;
  AnimatedTransform worldToPrimitive;
}
