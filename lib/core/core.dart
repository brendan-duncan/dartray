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
library core;

import 'dart:async';
import 'dart:math' as Math;
import 'dart:typed_data';
import 'package:archive/archive.dart';
import 'package:image/image.dart' as Img;

part 'light/area_light.dart';
part 'light/light_sample.dart';
part 'light/light_sample_offsets.dart';
part 'light/shape_set.dart';
part 'light/visibility_tester.dart';
part 'primitive/aggregate.dart';
part 'primitive/geometric_primitive.dart';
part 'primitive/transformed_primitive.dart';
part 'reflection/anisotropic.dart';
part 'reflection/blinn.dart';
part 'reflection/brdf_remap.dart';
part 'reflection/brdf_to_btdf.dart';
part 'reflection/bsdf.dart';
part 'reflection/bsdf_sample.dart';
part 'reflection/bsdf_sample_offsets.dart';
part 'reflection/bssrdf.dart';
part 'reflection/bxdf.dart';
part 'reflection/fresnel.dart';
part 'reflection/fresnel_blend.dart';
part 'reflection/fresnel_conductor.dart';
part 'reflection/fresnel_dielectric.dart';
part 'reflection/fresnel_no_op.dart';
part 'reflection/irregular_isotropic_brdf.dart';
part 'reflection/lambertian.dart';
part 'reflection/microfacet.dart';
part 'reflection/microfacet_distribution.dart';
part 'reflection/oren_nayar.dart';
part 'reflection/regular_halfangle_brdf.dart';
part 'reflection/scaled_bxdf.dart';
part 'reflection/specular_reflection.dart';
part 'reflection/specular_transmission.dart';
part 'texture/constant_texture.dart';
part 'texture/cylindrical_mapping_2d.dart';
part 'texture/identity_mapping_3d.dart';
part 'texture/planar_mapping_2d.dart';
part 'texture/spherical_mapping_2d.dart';
part 'texture/texture_mapping_2d.dart';
part 'texture/texture_mapping_3d.dart';
part 'texture/uv_mapping_2d.dart';
part 'volume/aggregate_volume.dart';
part 'volume/density_region.dart';
part 'volume/volume.dart';
part 'volume/volume_region.dart';
part 'animated_transform.dart';
part 'bbox.dart';
part 'camera.dart';
part 'camera_sample.dart';
part 'common.dart';
part 'differential_geometry.dart';
part 'film.dart';
part 'filter.dart';
part 'light.dart';
part 'integrator.dart';
part 'intersection.dart';
part 'kdtree.dart';
part 'log.dart';
part 'material.dart';
part 'matrix4x4.dart';
part 'mipmap.dart';
part 'montecarlo.dart';
part 'normal.dart';
part 'octree.dart';
part 'output_image.dart';
part 'param_set.dart';
part 'pixel_sampler.dart';
part 'plugin.dart';
part 'point.dart';
part 'primitive.dart';
part 'projective_camera.dart';
part 'quaternion.dart';
part 'ray.dart';
part 'ray_differential.dart';
part 'renderer.dart';
part 'render_overrides.dart';
part 'resource_manager.dart';
part 'rgb_color.dart';
part 'rng.dart';
part 'sample.dart';
part 'sampled_spectrum.dart';
part 'sampler.dart';
part 'scene.dart';
part 'shape.dart';
part 'spectrum.dart';
part 'spectrum_image.dart';
part 'spherical_harmonics.dart';
part 'stats.dart';
part 'surface_integrator.dart';
part 'surface_point.dart';
part 'texture.dart';
part 'texture_params.dart';
part 'transform.dart';
part 'vector.dart';
part 'volume_integrator.dart';
part 'xyz_color.dart';
