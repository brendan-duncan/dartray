
const String SCENE = """
SurfaceIntegrator "path" "integer maxdepth" [2]
Sampler "lowdiscrepancy" "integer pixelsamples" [2]

LookAt 0 0 -35 0 0 0 0 1 0 
Camera "perspective" "float fov" [35]

WorldBegin
############################################################################
# Light Source Definitions
############################################################################

AttributeBegin
AreaLightSource "area" "color L" [36 36 36 ] "integer nsamples" [1]
Translate 0 9.9 0
Rotate 90 1 0 0
Shape "disk" "float radius" [3]
AttributeEnd

############################################################################
# Wall Definitions
############################################################################
AttributeBegin
  Material "matte" "color Kd" [0.75 0.75 0.75 ] 
    Shape "trianglemesh"  "integer indices" [0 1 2 0 2 3 ] "point P" [ 10 -10 -10 -10 -10 -10 -10 -10  10  10 -10  10 ] 
    Shape "trianglemesh"  "integer indices" [0 1 2 0 2 3 ] "point P" [ 10  10 -10  10  10  10 -10  10  10 -10  10 -10 ] 
    Shape "trianglemesh"  "integer indices" [0 1 2 0 2 3 ] "point P" [ 10 -10  10 -10 -10  10 -10  10  10  10  10  10 ] 
  Material "matte" "color Kd" [0.48 0.1125 0.075]
    Shape "trianglemesh"  "integer indices" [0 1 2 0 2 3 ] "point P" [-10 -10  10 -10 -10 -10 -10  10 -10 -10  10  10 ] 
  Material "matte" "color Kd" [0.1125 0.375 0.1125]
    Shape "trianglemesh"  "integer indices" [0 1 2 0 2 3 ] "point P" [ 10 -10 -10  10 -10  10  10  10  10  10  10 -10 ] 
AttributeEnd

############################################################################
# Short Box Definition
############################################################################
AttributeBegin
  Translate 4 -7 4
  Scale 0.3 0.4 0.3
  Rotate 30 0 1 0
  #Material "glass" "float index" [1.5 ] "color Kr" [0.1 0.8 0.8 ] "color Kt" [0.1 0.8 0.8 ]
  Material "shinymetal" "float roughness" [10 ] "color Kr" [0.4 0.4 0 ] "color Ks" [0.6 0.6 0 ]
  Shape "trianglemesh"  "integer indices" [0 2 1 0 3 2 ] "point P" [ 10 -10 -10 -10 -10 -10 -10 -10  10  10 -10  10 ] 
  Shape "trianglemesh"  "integer indices" [0 2 1 0 3 2 ] "point P" [ 10  10 -10  10  10  10 -10  10  10 -10  10 -10 ] 
  Shape "trianglemesh"  "integer indices" [0 2 1 0 3 2 ] "point P" [ 10 -10  10 -10 -10  10 -10  10  10  10  10  10 ] 
  Shape "trianglemesh"  "integer indices" [0 2 1 0 3 2 ] "point P" [-10 -10  10 -10 -10 -10 -10  10 -10 -10  10  10 ] 
  Shape "trianglemesh"  "integer indices" [0 2 1 0 3 2 ] "point P" [ 10 -10 -10  10 -10  10  10  10  10  10  10 -10 ] 
  Shape "trianglemesh"  "integer indices" [0 1 2 0 2 3 ] "point P" [ 10 -10 -10 -10 -10 -10 -10  10 -10  10  10 -10 ] 
AttributeEnd

############################################################################
# Glass Sphere Definition
############################################################################
AttributeBegin
  #Material "shinymetal" "float roughness" [10 ] "color Kr" [0.4 0.4 0 ] "color Ks" [0.6 0.6 0 ]
  Material "glass" "float index" [1.5 ] "color Kr" [0.1 0.8 0.8 ] "color Kt" [0.1 0.8 0.8 ] 
  Translate -4 -4 0
    Shape "sphere" "float radius" 3  
AttributeEnd

WorldEnd
""";
