part of volume_regions;

class ExponentialDensityRegion extends DensityRegion {
  static ExponentialDensityRegion Create(Transform volume2world, ParamSet params) {
    // Initialize common volume region parameters
    Spectrum sigma_a = params.findOneSpectrum("sigma_a", new Spectrum(0.0));
    Spectrum sigma_s = params.findOneSpectrum("sigma_s", new Spectrum(0.0));
    double g = params.findOneFloat("g", 0.0);
    Spectrum Le = params.findOneSpectrum("Le", new Spectrum(0.0));
    Point p0 = params.findOnePoint("p0", new Point(0.0, 0.0, 0.0));
    Point p1 = params.findOnePoint("p1", new Point(1.0, 1.0, 1.0));
    double a = params.findOneFloat("a", 1.0);
    double b = params.findOneFloat("b", 1.0);
    Vector up = params.findOneVector("updir", new Vector(0.0, 1.0, 0.0));

    return new ExponentialDensityRegion(sigma_a, sigma_s, g, Le,
                                        new BBox(p0, p1),
                                        volume2world, a, b, up);
  }

  ExponentialDensityRegion(Spectrum sa, Spectrum ss,
                     double gg, Spectrum emit, this.extent,
                     Transform v2w, this.a, this.b,
                     Vector up) :
    super(sa, ss, gg, emit, v2w) {
      upDir = Vector.Normalize(up);
  }

  BBox worldBound() {
    return Transform.Inverse(worldToVolume).transformBBox(extent);
  }

  bool intersectP(Ray r, List<double> t0, List<double> t1) {
    Ray ray = worldToVolume.transformRay(r);
    return extent.intersectP(ray, t0, t1);
  }

  double density(Point Pobj) {
    if (!extent.inside(Pobj)) {
      return 0.0;
    }

    double height = Vector.Dot(Pobj - extent.pMin, upDir);

    return a * exp(-b * height);
  }

  BBox extent;
  double a;
  double b;
  Vector upDir;
}
