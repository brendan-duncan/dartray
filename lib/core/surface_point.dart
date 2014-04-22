part of core;

class SurfacePoint {
  SurfacePoint(Point p, Normal n, this.area, this.rayEpsilon) :
    this.p = new Point.from(p),
    this.n = new Normal.from(n);

  final Point p;
  final Normal n;
  final double area;
  final double rayEpsilon;
}
