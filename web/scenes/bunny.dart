
const String SCENE = """
Sampler "lowdiscrepancy" "integer pixelsamples" [4] 
SurfaceIntegrator "path"
#SurfaceIntegrator "ambientocclusion"
LookAt 0 .2 .2    -.02 .1 0  0 1 0
Camera "perspective" "float fov" [60]

WorldBegin

#AttributeBegin
#CoordSysTransform "camera"
#LightSource "point" "color I" [ 150 150 150 ] "point from" [5 10 5]
#AttributeEnd

AttributeBegin
  AreaLightSource "area" "color L" [100 100 100] "integer nsamples" [3]
  Translate 0 10 0
  Rotate 45 0 1 0
  Shape "trianglemesh" "integer indices" [ 0 1 2 2 3 0 ]
      "point P" [-.6 0 -.6   .6 0 -.6   .6 0 .6   -.6 0 .6 ]
AttributeEnd

Material "matte" "color Kd" [.4 .82 .4]
Shape "trianglemesh" "point P" [ -1 0 -1 1 0 -1 1 0 1 -1 0 1 ]
  "integer indices" [ 0 1 2 2 3 0]

Material "matte" "color Kd" [.6 .42 .4]

Include "geometry/bunny.pbrt"

WorldEnd
""";
