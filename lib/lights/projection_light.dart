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
part of lights;

class ProjectionLight extends Light {
  ProjectionLight(Transform light2world, this.intensity, String texname,
                  double fov)
      : super(light2world) {
    lightPos = lightToWorld.transformPoint(new Point(0.0, 0.0, 0.0));

    if (texname.isNotEmpty) {
      Completer completer = new Completer();
      ResourceManager.RequestImage(texname, completer.future)
        .then((SpectrumImage img) {
          String name = MIPMap.GetTextureName(texname);

          if (ResourceManager.HasTexture(name)) {
            projectionMap = ResourceManager.GetTexture(name);
          } else {
            projectionMap = new MIPMap.texture(img, texname);
            ResourceManager.AddTexture(name, projectionMap);
          }

          // Initialize projection matrix using the aspect ratio of the loaded
          // map.
          double aspect = img.width / img.height;
          if (aspect > 1.0)  {
            screenX0 = -aspect;
            screenX1 = aspect;
            screenY0 = -1.0;
            screenY1 = 1.0;
          } else {
            screenX0 = -1.0;
            screenX1 = 1.0;
            screenY0 = -1.0 / aspect;
            screenY1 = 1.0 / aspect;
          }

          // Compute cosine of cone surrounding projection directions
          double opposite = Math.tan(Radians(fov) / 2.0);
          double tanDiag = opposite * Math.sqrt(1.0 + 1.0 / (aspect * aspect));
          cosTotalWidth = Math.cos(Math.atan(tanDiag));

          completer.complete();
        });
    }

    // Initialize projection matrix
    double aspect = 1.0;
    if (aspect > 1.0)  {
      screenX0 = -aspect;
      screenX1 = aspect;
      screenY0 = -1.0;
      screenY1 = 1.0;
    } else {
      screenX0 = -1.0;
      screenX1 = 1.0;
      screenY0 = -1.0 / aspect;
      screenY1 = 1.0 / aspect;
    }
    hither = 1.0e-3;
    yon = 1.0e30;
    lightProjection = Transform.Perspective(fov, hither, yon);

    // Compute cosine of cone surrounding projection directions
    double opposite = Math.tan(Radians(fov) / 2.0);
    double tanDiag = opposite * Math.sqrt(1.0 + 1.0 / (aspect * aspect));
    cosTotalWidth = Math.cos(Math.atan(tanDiag));
  }

  Spectrum sampleLAtPoint(Point p, double pEpsilon, LightSample ls,
                          double time, Vector wi, List<double> pdf,
                          VisibilityTester visibility) {
    wi.copy(Vector.Normalize(lightPos - p));
    pdf[0] = 1.0;
    visibility.setSegment(p, pEpsilon, lightPos, 0.0, time);
    return intensity * projection(-wi) / Vector.DistanceSquared(lightPos, p);
  }

  bool isDeltaLight() {
    return true;
  }

  Spectrum projection(Vector w) {
    Vector wl = worldToLight.transformVector(w);
    // Discard directions behind projection light
    if (wl.z < hither) {
      return new Spectrum(0.0);
    }

    // Project point onto projection plane and compute light
    Point Pl = lightProjection.transformPoint(new Point.from(wl));
    if (Pl.x < screenX0 || Pl.x > screenX1 ||
        Pl.y < screenY0 || Pl.y > screenY1) {
      return new Spectrum(0.0);
    }

    if (projectionMap == null) {
      return new Spectrum(1.0);
    }

    double s = (Pl.x - screenX0) / (screenX1 - screenX0);
    double t = (Pl.y - screenY0) / (screenY1 - screenY0);

    return new Spectrum.from(projectionMap.lookup(s, t),
                             Spectrum.SPECTRUM_ILLUMINANT);
  }

  Spectrum power(Scene scene) {
    return (projectionMap != null ?
            new Spectrum.from(projectionMap.lookup(0.5, 0.5, 0.5),
                              Spectrum.SPECTRUM_ILLUMINANT) :
            new Spectrum(1.0)) *
            intensity * (2.0 * Math.PI * (1.0 - cosTotalWidth));
  }

  Spectrum sampleL(Scene scene, LightSample ls, double u1, double u2,
                   double time, Ray ray, Normal Ns, List<double> pdf) {
    Vector v = UniformSampleCone(ls.uPos[0], ls.uPos[1], cosTotalWidth);
    ray.set(lightPos, lightToWorld.transformVector(v), 0.0, INFINITY, time);
    Ns.copy(ray.direction);
    pdf[0] = UniformConePdf(cosTotalWidth);
    return intensity * projection(ray.direction);
  }

  double pdf(Point p, Vector w) {
    return 0.0;
  }

  static ProjectionLight Create(Transform light2world, ParamSet paramSet) {
    Spectrum I = paramSet.findOneSpectrum('I', new Spectrum(1.0));
    Spectrum sc = paramSet.findOneSpectrum('scale', new Spectrum(1.0));
    double fov = paramSet.findOneFloat('fov', 45.0);
    String texname = paramSet.findOneFilename('mapname', '');
    return new ProjectionLight(light2world, I * sc, texname, fov);
  }

  MIPMap projectionMap;
  Point lightPos;
  Spectrum intensity;
  Transform lightProjection;
  double hither;
  double yon;
  double screenX0;
  double screenX1;
  double screenY0;
  double screenY1;
  double cosTotalWidth;
}
