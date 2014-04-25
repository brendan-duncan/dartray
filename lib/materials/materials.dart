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
library materials;

import 'dart:async';
import 'dart:typed_data';
import 'dart:math' hide Point;
import '../core/core.dart';

import 'package:image/image.dart';

part 'glass_material.dart';
part 'kd_subsurface_material.dart';
part 'matte_material.dart';
part 'measured_material.dart';
part 'metal_material.dart';
part 'mirror_material.dart';
part 'mix_material.dart';
part 'plastic_material.dart';
part 'shiny_metal_material.dart';
part 'substrate_material.dart';
part 'subsurface_material.dart';
part 'translucent_material.dart';
part 'uber_material.dart';
