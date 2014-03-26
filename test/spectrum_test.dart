part of dartray_test;

void spectrum_test() {
  group('RGBColor', () {
    test('constructor', () {
      RGBColor c1 = new RGBColor();
      expect(c1.c[0], equals(0.0));
      expect(c1.c[1], equals(0.0));
      expect(c1.c[2], equals(0.0));
      expect(c1.isBlack(), equals(true));

      RGBColor c2 = new RGBColor(3.0);
      expect(c2.c[0], equals(3.0));
      expect(c2.c[1], equals(3.0));
      expect(c2.c[2], equals(3.0));
      expect(c2.isBlack(), equals(false));

      RGBColor c3 = new RGBColor.rgb(0.0, 1.0, 2.0);
      expect(c3.c[0], equals(0.0));
      expect(c3.c[1], equals(1.0));
      expect(c3.c[2], equals(2.0));
      expect(c3.isBlack(), equals(false));

      RGBColor c4 = new RGBColor.from(c3);
      expect(c4.c[0], equals(0.0));
      expect(c4.c[1], equals(1.0));
      expect(c4.c[2], equals(2.0));
    });
  });
}
