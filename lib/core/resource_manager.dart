/****************************************************************************
 *  Copyright (C) 2014 by Brendan Duncan.                                   *
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
 *   pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.*
 ****************************************************************************/
part of core;

abstract class ResourceManager {
  static ResourceManager global;

  List<String> includePaths = [];

  ResourceManager() {
    global = this;
  }

  static Future<List<int>> RequrestFile(String path) {
    return global.requrestFile(path);
  }

  static List<int> GetFile(String path) {
    return global.getFile(path);
  }

  void addIncludePath(path) {
    includePaths.add(path);
  }

  Future<List<int>> requrestFile(String path);

  List<int> getFile(String path) {
    if (!resources.containsKey(path)) {
      return null;
    }
    return resources[path];
  }

  Map<String, List<int>> resources = {};
}
