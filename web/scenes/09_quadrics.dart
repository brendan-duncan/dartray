const String SCENE = """
#Sampler "bestcandidate" "integer pixelsamples" [64]
PixelFilter "mitchell"

Scale -1 1 1
LookAt 0 3 8  0 .8 0   0 1 0
#Scale -1 1 1
Camera "perspective" "float fov" [21]

SurfaceIntegrator "directlighting"

WorldBegin

AttributeBegin
  CoordSysTransform "camera"
  AreaLightSource "area" "color L" [8 8 8  ] "integer nsamples" [4]
  Translate 0 2 -10
  Material "matte" "color Kd" [ 0 0 0 ]
  Shape "disk" "float radius" [3] 
AttributeEnd

AttributeBegin
  AreaLightSource "area" "color L" [3.2 3.2 3.2] "integer nsamples" [4]
  Translate 0 10 0
  Rotate 90 1 0 0 
  Material "matte" "color Kd" [ 0 0 0 ]
  Shape "disk" "float radius" [20] 
AttributeEnd

AttributeBegin
  Material "matte" "color Kd" [.8 .8 .8 ]
  Shape "trianglemesh" "integer indices" [ 0 1 2 2 0 3 ]
    "point P" [-10 0 -10   10 0 -10   10 0 10   -10 0 10 ]
  Shape "trianglemesh" "integer indices" [ 0 1 2 2 0 3 ]
    "point P" [-10 0 -10   10 0 -10   10 9 -10   -10 9 -10 ]
  Shape "trianglemesh" "integer indices" [ 0 1 2 2 0 3 ]
    "point P" [-10 0 -10   -10 0 10   -10 9 10   -10 9 -10 ]
AttributeEnd

Texture "g" "color" "imagemap" "string filename" "textures/grid.exr"
Texture "grid" "color" "scale" "texture tex1" "g" "color tex2" [.75 .75 .75]

Material "uber" "texture Kd" "grid" "color Ks" [.1 .1 .1]
  "color opacity" [.8 .8 .8]

AttributeBegin
Translate -1.75 0 -.4
Scale .7 1.8 .7
Rotate 60 0 1 0
Rotate -90 1 0 0
Shape "paraboloid"
AttributeEnd

AttributeBegin
Translate 1.75 0 -.4
Scale .7 1.8 .7
Rotate 54 0 1 0
Rotate -90 1 0 0
Shape "cone"
AttributeEnd

AttributeBegin
Translate -.15 1.35 -.4
Scale .5 1.7 .5
Rotate 150 0 1 0
Rotate 90 1 0 0 
Shape "hyperboloid" "point p1" [ 1 0 0 ] "point p2" [ .8 1 .8 ] 
AttributeEnd

WorldEnd
""";
