const String SCENE = """
Sampler "stratified" 
  "integer xsamples" [5] "integer ysamples" [5]
  "bool jitter" ["false"]

PixelFilter "mitchell"

#Accelerator "none"
Accelerator "grid"

LookAt -4 3 -11.3  -4 3 0   0 1 0
Camera "perspective"
  "float fov" [60]

#SurfaceIntegrator "directlighting" "integer maxdepth" [6]
SurfaceIntegrator "whitted" "integer maxdepth" [6]

WorldBegin

# ----------
# Shapes
# ----------
Texture "floor" "spectrum" "checkerboard" "float uscale" 16 "float vscale" 24 "string aamode" "closedform"
  "color tex2" [.82 .16 .07] "color tex1" [.88 .89 .0]
AttributeBegin # floor
  Material "matte" "texture Kd" "floor"
  Shape "trianglemesh" "integer indices" [ 0 1 2 2 0 3 ]
    "point P" [-8 0 -8   8 0 -8   8 0 16   -8 0 16 ]
    "float uv" [0 0   1 0   1 1   0 1]
AttributeEnd

AttributeBegin # white sphere behind
  Material "uber" "color Kd" [.8 .8 .8 ] "color Ks" [.8 .8 .8 ]  "color Kr" [.8 .8 .8 ]  "color Kt" [0 0 0]  "float roughness" [0.005] 
    Translate -2.4 2.3 -5.1
  Shape "sphere" "float radius" [ 1.15]
AttributeEnd

AttributeBegin # glass sphere in front
  Material "uber" 
  "color Kd" [.2 .2 .2 ] 
  "color Ks" [1 1 1 ]  
  "color Kr" [1 1 1 ]  
  "color Kt" [.9 .9 .9]  
  "float roughness" [0] "float index" [1.05] 
    Translate -4.2 3.1-6.6
    Shape "sphere" "float radius" [ 1.3]
AttributeEnd



# ----------
# Lights
# ----------
AttributeBegin # the sun
  LightSource "distant" "blackbody L" [3200 1] 
  "point from" [-2 15 -10] "point to" [0 0 0] "integer nsamples" [1]
AttributeEnd


AttributeBegin # the sky
  LightSource "infinite" "blackbody L" [20000 .8] "integer nsamples" [1]
AttributeEnd


WorldEnd
""";
