/****************************************************************************
 *  Copyright (C) 2014 by authors (see AUTHORS)                             *
 *                                                                          *
 *  This file is part of DartRay.                                           *
 *                                                                          *
 *  Licensed under the Apache License, Version 2.0 (the "License");         *
 *  you may not use this file except in compliance with the License.        *
 *  You may obtain a copy of the License at                                 *
 *                                                                          *
 *  http://www.apache.org/licenses/LICENSE-2.0                              *
 *                                                                          *
 *  Unless required by applicable law or agreed to in writing, software     *
 *  distributed under the License is distributed on an "AS IS" BASIS,       *
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 *  See the License for the specific language governing permissions and     *
 *  limitations under the License.                                          *
 *                                                                          *
 *   This project is based on PBRT v2 ; see http://www.pbrt.org             *
 ****************************************************************************/
library accelerators;

import 'dart:math' as Math;
import 'dart:typed_data';
import '../core/core.dart';

part 'bvh_accel.dart';
part 'grid_accel.dart';
part 'kdtree_accel.dart';
part 'naive_accel.dart';
