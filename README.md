# DartRay

##Overview

DartRay is a free, open-source physically based renderer, providing a variety
of rendering, lighting and material algorithms.

DartRay port of the PBRT renderer (www.pbrt.org).

##Caveats
DartRay provides a test-bed for various rendering algorithms. 
The primary (and perhaps only) objectives of DartRay are experimentation, 
learning, and most important, fun. While it is perfectly capable of producing 
incredible images, it has a lot of things going against it to be thought of as 
anything remotely close to production worthy. 

For one, a production renderer should not be written in an interpreted 
language, even one as great as Dart. Dart, being a web-centric language, does
have some limitations that make scalability more difficult.

Dart supports multiple threads using a mechanism called isolates. An isolate
is more analogous to multi-processing than multi-threading. Each isolate 
has its own memory space, and cannot share any memory with any other
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
a great environment for experimentation and learning about rendering techniques.
It provides a reasonably fast and flexible, cross-platform environment for
experimentation and learning. And, most importantly, it's fun. These are the 
goals of DartRay.

##License
Copyright (C) 2014 Brendan Duncan

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

----

DartRay is a derivative of the PBRT-V2 Rendering System by Greg Humphreys and
Matt Pharr.

Copyright (c) 1998-2012, Matt Pharr and Greg Humphreys.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list 
of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this 
list of conditions and the following disclaimer in the documentation and/or 
other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
