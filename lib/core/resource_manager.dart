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

  static Future<List<int>> RequestFile(String path, [Future future]) {
    return global.requestFile(path, future);
  }

  static Future<String> RequestScene(String path, [Future future]) {
    return global.requestScene(path, future);
  }

  static Future<SpectrumImage> RequestImage(String path, [Future future]) {
    return global.requestImage(path, future);
  }

  static List<int> GetFile(String path) {
    return global.getFile(path);
  }

  void addIncludePath(path) {
    includePaths.add(path);
  }

  bool hasFile(String path) => resources.containsKey(path);

  Future<List<int>> requestFile(String path, [Future future]);

  Future<String> requestScene(String path, [Future future]) {
    if (resources.containsKey(path)) {
      if (resources[path] is Future) {
        return resources[path];
      }

      Completer<List<int>> c = new Completer<List<int>>();
      c.complete(resources[path]);
      return c.future;
    }

    Completer<String> c = new Completer<String>();
    resources[path] = c.future;

    requestFile(path, future).then((bytes) {
      print('SCENE $path LOADED');
      String s = new String.fromCharCodes(bytes);
      resources[path] = s;
      c.complete(s);
    });
    return c.future;
  }

  Future<SpectrumImage> requestImage(String path, [Future future]) {
    if (resources.containsKey(path)) {
      if (resources[path] is Future) {
        return resources[path];
      }

      Completer<List<int>> c = new Completer<List<int>>();
      c.complete(resources[path]);
      return c.future;
    }

    Completer<SpectrumImage> c = new Completer<SpectrumImage>();
    resources[path] = c.future;

    requestFile(path, future).then((bytes) {
      print('IMAGE $path LOADED');
      SpectrumImage img = new SpectrumImage(1, 1);
      resources[path] = img;
      c.complete(img);
    });

    return c.future;
  }

  Future waitUntilReady() {
    Completer c = new Completer();
    Future.wait(futures).then((List responses) {
      futures.clear();
      c.complete();
    });
    return c.future;
  }

  getFile(String path) {
    if (!resources.containsKey(path)) {
      return null;
    }
    return resources[path];
  }

  List<Future> futures = [];
  Map<String, dynamic> resources = {};
}
