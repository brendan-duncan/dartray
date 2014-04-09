const String SCENE = """
Sampler "stratified" 
  "integer xsamples" [2] "integer ysamples" [2]
  "bool jitter" ["false"]

PixelFilter "mitchell"

Accelerator "grid"

LookAt 2 4 10   0 0 0   0 1 0

Camera "perspective"
  "float fov" [45]

SurfaceIntegrator "directlighting"
  "integer maxdepth" [5]

WorldBegin

# ----------
# Shapes
# ----------
AttributeBegin # large white disk
  Rotate 90 1 0 0
  Material "matte" "color Kd" [1 1 1]
  Shape "disk" "float height" [1] "float radius" [5]  
AttributeEnd

AttributeBegin # red box
  Translate 2.1 0 0
  Material "matte" "color Kd" [1 0 0]
  Shape "trianglemesh"
    "integer indices"  [
      0 1 2  0 2 3 # back
      0 5 1 0 4 5 # up
      1 5 6 1 6 2 #right
      4 6 5 4 7 6 # front
      0 7 3 0 4 7 # left
      2 6 7 2 7 3 # down
       ] 
    "point P" [-1 1 1  1 1 1  1 -1 1 -1 -1 1  -1 1 -1  1 1 -1  1 -1 -1 -1 -1 -1]
AttributeEnd

AttributeBegin # yellow sphere
  Translate 0 0 -2.1
  Material "matte" "color Kd" [0.9 0.8 0]
  Shape "sphere"  
AttributeEnd

AttributeBegin # blue box
  Translate -2.1 0 0
  Material "matte" "color Kd" [0 0 1]
  Shape "trianglemesh"
    "integer indices"  [
      0 1 2  0 2 3 # back
      0 5 1 0 4 5 # up
      1 5 6 1 6 2 #right
      4 6 5 4 7 6 # front
      0 7 3 0 4 7 # left
      2 6 7 2 7 3 # down
       ] 
    "point P" [-1 1 1  1 1 1  1 -1 1 -1 -1 1  -1 1 -1  1 1 -1  1 -1 -1 -1 -1 -1] 
AttributeEnd

# ----------
# Lights
# ----------
AttributeBegin
  Translate 0 20 -20
  LightSource "point" "color I" [1200 1200 1200] 
AttributeEnd

AttributeBegin
  Translate 20 20 -20
  LightSource "point" "color I" [800 800 800]
AttributeEnd

AttributeBegin
  Translate -20 -2 1
  LightSource "point" "color I" [160 160 240]  
AttributeEnd

AttributeBegin
  Translate 20 -2 1
  LightSource "point" "color I" [260 260 320] 
AttributeEnd

AttributeBegin
  Translate 0 2 4
  LightSource "point" "color I" [20 20 20] 
AttributeEnd


WorldEnd
""";
