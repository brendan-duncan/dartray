# DartRay

##Overview

DartRay is an extendable physically based ray tracer providing a number of 
state-of-the-art algorithms for rendering, lighting and shading. It is written 
in the [Dart](www.dartlang.org) programming language, with an emphasis on experimentation
and learning.

DartRay is a port of the [PBRT](http://www.pbrt.org) renderer by Matt Pharr and 
Greg Humphreys. I highly recommend the book on PBRT, [Physically Based Rendering: From Theory To Implementation](http://www.amazon.com/gp/product/0123750792?ie=UTF8&tag=pharr-20&linkCode=as2&camp=1789&creative=390957&creativeASIN=0123750792).

##[Documentation](https://github.com/brendan-duncan/dartray/wiki)

##[Test DartRay](http://brendan-duncan.github.io/dartray/web_dartray/web_dartray.html)<br>
*_Note that DartRay was written for the Dart VM and Javascript will be significantly slower._*


##Test Renders
![Dipole Subsurface Scattering](https://cloud.githubusercontent.com/assets/3642099/2797718/917d9174-cc39-11e3-881d-3b8f16d10be1.png)
![Path Tracing surface integrator](https://cloud.githubusercontent.com/assets/3642099/2703164/cf037684-c447-11e3-81d7-e2f2c6520fa1.jpg)
![IGI (Instant Global Illumination) surface integrator, environment camera](https://cloud.githubusercontent.com/assets/3642099/2714447/0e5122c8-c4f7-11e3-8436-a7c5b9011ab7.jpg)
![Ambient Occlusion surface integrator](https://cloud.githubusercontent.com/assets/3642099/2650559/0416dcb2-bf75-11e3-8c3e-fdb836e8146b.jpg)
![smoke](https://cloud.githubusercontent.com/assets/3642099/2797791/db44e606-cc3c-11e3-8617-beb25d29dfbc.png)
![yeahright](https://cloud.githubusercontent.com/assets/3642099/2804809/abddafde-ccae-11e3-8a83-6089a06e5e1b.png)

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
