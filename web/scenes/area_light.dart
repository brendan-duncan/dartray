
const String SCENE = """
SurfaceIntegrator "path"
Sampler "lowdiscrepancy" "integer pixelsamples" [16] 
Accelerator "bvh"
#Accelerator "grid"

LookAt 4 4 10   0 2 0  0 1 0
Camera "perspective" "float fov" [30]

WorldBegin
  AttributeBegin
    AreaLightSource "area" "color L" [100 100 100] "integer nsamples" [2]
    Translate 1 4.5 1
    Rotate 45 0 1 0
    Shape "trianglemesh" "integer indices" [ 0 1 2 2 3 0 ]
        "point P" [-.6 0 -.6   .6 0 -.6   .6 0 .6   -.6 0 .6 ]
  AttributeEnd

  AttributeBegin
    Material "matte" "color Kd" [.5 .5 .5]
    Translate 0 1 0

    Shape "trianglemesh" 
    #Shape "loopsubdiv" "integer nlevels" [2]
      "integer indices" [
          0 1 2
          0 2 3
          4 5 6
          4 6 7
          0 4 7
          0 7 1
          1 7 6
          1 6 2
          2 6 5
          2 5 3
          4 0 3
          4 3 5] 
      "point P" [
          1.000000 1.000000 -1.000000
          1.000000 -1.000000 -1.000000
          -1.000000 -1.000000 -1.000000
          -1.000000 1.000000 -1.000000
          1.000000 0.999999 1.000000
          -1.000000 1.000000 1.000000
          -1.000000 -1.000000 1.000000
          0.999999 -1.000001 1.000000]
  AttributeEnd

  AttributeBegin
    #Material "matte" "color Kd" [.8 .8 .8]
    Material "uber" "color Kd" [0 0 1] "color opacity" [.5 .5 .5] "color Kr" [0 0 1] "color Kt" [0 0 0.5]
    Shape "trianglemesh" 
        "point P" [ -100 0 -100   100 0 -100   100 0 100   -100 0 100 ]
        "float uv" [ 0 0 1 0 1 1 0 1 ]
        "integer indices" [ 0 1 2 2 3 0]
  AttributeEnd
WorldEnd
""";
