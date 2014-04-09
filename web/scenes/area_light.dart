
const String SCENE = """
SurfaceIntegrator "path"
#Sampler "lowdiscrepancy" "integer pixelsamples" [4] 
#Sampler "stratified" "integer pixelsamples" [1]
Sampler "random" "integer pixelsamples" [1]
#Sampler "adaptive" "integer minSamples" [2] "integer maxSamples" [4]
#Sampler "halton"
#Sampler "bestcandidate"

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
    Material "matte" "color Kd" [.6 .4 .4]
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
    Material "matte" "color Kd" [.4 .6 .4]
    Shape "trianglemesh" 
        "point P" [ -100 0 -100   100 0 -100   100 0 100   -100 0 100 ]
        "float uv" [ 0 0 1 0 1 1 0 1 ]
        "integer indices" [ 0 1 2 2 3 0]
  AttributeEnd
WorldEnd
""";
