# image

##Overview

DartRay is a free, open-source physically based renderer, providing a variety
of rendering algorithms, lighting and material models.

DartRay port of the PBRT renderer (www.pbrt.org).

DartRay provides a test-bed for various rendering algorithms. 
The primary (and perhaps only) objective of DartRay is experimentation and 
learning. While it is perfectly capable of producing incredible images, it has 
a lot of things going against it to be remotely thought of as anything close to 
production worthy. 

For one, a production renderer should not be written in an interpreted 
language, even one as great as Dart. Dart, being a web-centric language, does
have some limitations that keep make scalability more difficult. 

Dart supports multiple threads using a mechanism called isolates. An isolate
is more anagolous to multi-processing than multi-threading. Each isolate 
(thread) has its own memory space, and cannot share any memory with any other
isolate. Communication between isolates is done through messaging, where the
data packet of the messages is copied as it's transfered to the other isolate.
What this means is that each isolate needs to load its own full copy of the 
scene and all of its resources. This will greatly restrict the number of threads
you can execute on large scenes, as well as increase the startup time for each 
isolate as they will each have to load all resources and preprocess the scene.  

Dart also lacks high-level GPGPU access such as OpenCL or Cuda.
Server side (dart:io) applications will not have any access to OpenGL, and
web applications wil have access to WebGL, but this diversity will make a 
uniform GPGPU acceleration impossible for both environments.

There are of course many other limitions. But on the positive side, Dart is
a great environment for experimentation and learing about rendering techniques.
It provides a reasonably fast and flexible, cross-platform environment for
experimentation and learning. And, most importantly, it's fun. These are the 
goals of DartRay.

