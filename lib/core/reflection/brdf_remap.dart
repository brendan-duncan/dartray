part of core;

Point BRDFRemap(Vector wo, Vector wi) {
  double cosi = Vector.CosTheta(wi);
  double coso = Vector.CosTheta(wo);
  double sini = Vector.SinTheta(wi);
  double sino = Vector.SinTheta(wo);
  double phii = Vector.SphericalPhi(wi);
  double phio = Vector.SphericalPhi(wo);
  double dphi = phii - phio;

  if (dphi < 0.0) {
    dphi += 2.0 * Math.PI;
  }

  if (dphi > 2.0 * Math.PI) {
    dphi -= 2.0 * Math.PI;
  }

  if (dphi > Math.PI) {
    dphi = 2.0 * Math.PI - dphi;
  }

  return new Point(sini * sino, dphi / Math.PI, cosi * coso);
}
