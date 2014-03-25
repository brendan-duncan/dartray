const String SCENE = """
Sampler "stratified" 
  "integer xsamples" [2] "integer ysamples" [2]
  "bool jitter" ["false"]

PixelFilter "mitchell"

Accelerator "grid"

LookAt 2 4 -6   0 0 0   0 1 0
Camera "perspective"
  "float fov" [50]

SurfaceIntegrator "directlighting"

WorldBegin

# ----------
# Shapes
# ----------
AttributeBegin # black sphere
  Translate 0 0 0
  Material "matte" "color Kd" [0.1 0.1 0.1]
  Shape "sphere" "float radius" [0.2]  
AttributeEnd
AttributeBegin # magenta disk
  Material "matte" "color Kd" [1 0 1]
  Shape "disk" "float height" [0] "float radius" [0.5]  
AttributeEnd
AttributeBegin # cyan disk
  Rotate 90 1 0 0
  Material "matte" "color Kd" [0 1 1]
  Shape "disk" "float height" [0] "float radius" [0.5]  
AttributeEnd
AttributeBegin # yellow disk
  Rotate 90 0 1 0
  Material "matte" "color Kd" [1 1 0]
  Shape "disk" "float height" [0] "float radius" [0.5]  
AttributeEnd
AttributeBegin # large white disk (floor)
  Rotate 90 1 0 0
  Material "matte" "color Kd" [1 1 1]
  Shape "disk" "float height" [1] "float radius" [50]  
AttributeEnd


AttributeBegin # red sphere
  Translate 0 0 2.1
  Material "matte" "color Kd" [1 0 0]
  Shape "sphere"  
AttributeEnd
AttributeBegin # green paraboloid
  Translate 2.1 1 0
  Rotate 90 1 0 0
  Scale 1 1 2
  Material "matte" "color Kd" [0 1 0]
  Shape "paraboloid"  
AttributeEnd

AttributeBegin # blue hyberboloid
  Translate -2.1 1 0
  Rotate 90 1 0 0
  Scale 1 1 2
  Material "matte" "color Kd" [0.1 0.1 1]
  Shape "hyperboloid"   "point p1" [ 1 0 0 ] "point p2" [ .8 1 1 ] 
AttributeEnd

AttributeBegin # yellow cone
  Translate 0 -1 -2.1
  Rotate -90 1 0 0
  Scale 1 1 2
  Material "matte" "color Kd" [0.9 0.7 0]
  Shape "cone" 
AttributeEnd

# ----------
# Lights
# ----------
AttributeBegin # sky sphere
  AreaLightSource "area" "color L" [.5 .8 1] "integer nsamples" [40] # 100 W Tungsten, 2600 K
  Translate 0 0 0
  Rotate -90 1 0 0
  ReverseOrientation
  Material "matte" "color Kd" [1 1 1]
  Shape "sphere"  "float radius" 50 "float zmin" -20
AttributeEnd
AttributeBegin # sun disk
  AreaLightSource "area" "color L" [2.04 1.578 1.144] "integer nsamples" [40] # 100 W Tungsten, 2600 K
  Rotate -20 1 0 0
  Translate 0 10 0
  Rotate 90 1 0 0
  Material "matte" "color Kd" [1 1 1]
  Shape "disk"  "float radius" 40
AttributeEnd

WorldEnd
""";
