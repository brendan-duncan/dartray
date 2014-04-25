/****************************************************************************
 * Copyright (C) 2014 by Brendan Duncan.                                    *
 *                                                                          *
 * This file is part of DartRay.                                            *
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License");          *
 * you may not use this file except in compliance with the License.         *
 * You may obtain a copy of the License at                                  *
 *                                                                          *
 * http://www.apache.org/licenses/LICENSE-2.0                               *
 *                                                                          *
 * Unless required by applicable law or agreed to in writing, software      *
 * distributed under the License is distributed on an "AS IS" BASIS,        *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 * See the License for the specific language governing permissions and      *
 * limitations under the License.                                           *
 *                                                                          *
 * This project is based on PBRT v2 ; see http://www.pbrt.org               *
 * pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.  *
 ****************************************************************************/
library surface_integrators;

import 'dart:async';
import 'dart:math' as Math;
import 'dart:typed_data';
import '../core/core.dart';
import '../renderers/renderers.dart';
import '../samplers/samplers.dart';

part 'ambient_occlusion_integrator.dart';
part 'diffuse_prt_integrator.dart';
part 'dipole_subsurface_integrator.dart';
part 'direct_lighting_integrator.dart';
part 'glossy_prt_integrator.dart';
part 'igi_integrator.dart';
part 'irradiance_cache_integrator.dart';
part 'path_integrator.dart';
part 'photon_map_integrator.dart';
part 'use_probes_integrator.dart';
part 'whitted_integrator.dart';
