library test2;

const String SCENE = """
#SurfaceIntegrator "path"
#Sampler "lowdiscrepancy" "integer pixelsamples" [2]

LookAt 4 4 10   0 2 0  0 1 0
Camera "perspective" "float fov" [30]

WorldBegin
  AttributeBegin
    Translate 0 10 10
    LightSource "point" "color I" [500 500 500]
  AttributeEnd

  AttributeBegin
    Material "metal"
    Translate 0 1 0
    Shape "sphere"
  AttributeEnd

  AttributeBegin
    Material "glass"
    Translate 0 3 0
    Shape "sphere"
  AttributeEnd

  AttributeBegin
    Material "shinymetal"
    Translate 2 1 0
    Shape "sphere"
  AttributeEnd

  AttributeBegin
    Material "translucent"
    Translate -2 1 0
    Shape "sphere"
  AttributeEnd
  
  AttributeBegin
    #Texture "tex" "color" "fbm"
    Texture "tex" "color" "checkerboard" "float uscale" [40]
            "float vscale" [40] "color tex1" [1 0 0] "color tex2" [0 0 1]
    Material "matte" "texture Kd" "tex"
    Shape "trianglemesh" 
        "point P" [ -100 0 -100   100 0 -100   100 0 100   -100 0 100 ]
        "float uv" [ 0 0 1 0 1 1 0 1 ]
        "integer indices" [ 0 1 2 2 3 0]
  AttributeEnd
WorldEnd
""";
