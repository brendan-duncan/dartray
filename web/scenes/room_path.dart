
const String SCENE = """
Scale -1 1 1
LookAt 18 3.5 2    8 -.5 20   0 1 0
Camera "perspective" "float fov" [50]

Sampler "lowdiscrepancy" "integer pixelsamples" [16]
#PixelFilter "box"

SurfaceIntegrator "path"

WorldBegin

# lights
LightSource "spot" "color I" [7000000 7000000 7000000] "point from" [70 230 -300]
    "point to" [10 10 10] "float coneangle" [1]
    "float conedeltaangle" [.01]

# walls
AttributeBegin
Material "matte" "color Kd" [.5 .5 .5]
#Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
#    "point P" [ 0 0 0    20 0 0    20 10 0    0 10 0 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ 0 0 20    20 0 20    20 10 20    0 10 20 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ 0 0 0    0 0 20    0 10 20    0 10 0 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ 20 0 0    20 0 20    20 10 20    20 10 0 ]

AttributeEnd

#ceiling
Material "matte" "color Kd" [.4 .4 .4]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ -3 10 -3  23 10 -3   23 10 5   -3 10 5 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ -3 10 15  23 10 15   23 10 23   -3 10 23 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ -3 10 -3  5 10 -3   5 10 23   -3 10 23 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ 15 10 -3  23 10 -3   23 10 23   15 10 23 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ -3 10.1 -3  23 10.1 -3   23 10.1 5   -3 10.1 5 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ -3 10.1 15  23 10.1 15   23 10.1 23   -3 10.1 23 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ -3 10.1 -3  5 10.1 -3   5 10.1 23   -3 10.1 23 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ 15 10.1 -3  23 10.1 -3   23 10.1 23   15 10.1 23 ]

# window bar things
Material "matte" "color Kd" [0 0 0 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ 7.4 10 0  7.6 10 0 7.6 10 20   7.4 10 20 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ 9.9 10 0  10.1 10 0 10.1 10 20   9.9 10 20 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ 12.4 10 0  12.6 10 0 12.6 10 20   12.4 10 20 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ 0 10 7.4  0 10 7.6 20 10 7.6   20 10 7.4 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ 0 10 9.9  0 10 10.1 20 10 10.1   20 10 9.9 ]
Shape "trianglemesh" "integer indices" [0 1 2 0 2 3 ]
    "point P" [ 0 10 12.4  0 10 12.6 20 10 12.6   20 10 12.4 ]

#floor
Material "matte" "color Kd" [.8 .8 .8]
Shape "trianglemesh" "integer indices" [0 1 2 0 3 2 ]
    "point P" [ 0 0 0  20 0 0   20 0 20   0 0 20 ]

WorldEnd
""";
