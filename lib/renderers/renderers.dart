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
library renderers;

import 'dart:async';
import 'dart:math' as Math;
import 'dart:typed_data';

import '../core/core.dart';
import '../pixel_samplers/pixel_samplers.dart';
import '../samplers/samplers.dart';
import '../shapes/shapes.dart';

part 'aggregate_test_renderer.dart';
part 'create_probes_renderer.dart';
part 'metropolis_renderer.dart';
part 'sampler_renderer.dart';
part 'surface_points_renderer.dart';
