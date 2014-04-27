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

  static Future<SpectrumImage> RequestImage(String path, [Future future]) {
    return global.requestImage(path, future);
  }

  static bool HasResource(String path) {
    return global.hasResource(path);
  }

  static GetResource(String path) {
    return global.getResource(path);
  }

  static void WriteFile(String path, data) {
    global.writeFile(path, data);
  }

  static bool HasTexture(String name) {
    return global.hasTexture(name);
  }

  static MIPMap GetTexture(String name) {
    return global.getTexture(name);
  }

  static void AddTexture(String name, MIPMap texture) {
    global.addTexture(name, texture);
  }

  /**
   * Add a path to look for resources in.
   */
  void addIncludePath(path) {
    includePaths.add(path);
  }

  /**
   * Has the resource with the given [path] been loaded?
   */
  bool hasResource(String path) => resources.containsKey(path);

  /**
   * Asynchronously load a file. This does not store the file as a resource,
   * so requesting the same file twice will load the file twice. Use the
   * resource based request methods instead, such as [requestTextFile] or
   * [requestImage].
   */
  Future<List<int>> loadFile(String path);

  /**
   * Request a binary file.
   *
   * Requesting the same file multiple times will only load the file once.
   * If the requester wants to process the contents of the file prior to
   * rendering starting, they must pass a [future] to the request. Once the
   * requester recieves the response from the request, it can then process
   * the data and, once finished, complete the [Completer] that the [future]
   * came from. When all requesters of data have completed, the render will
   * continue.
   */
  Future requestFile(String path, [Future future]) {
    if (future != null) {
      futures.add(future);
    }

    if (resources.containsKey(path)) {
      if (resources[path] is Future) {
        return resources[path];
      }

      Completer c = new Completer();
      c.complete(resources[path]);

      return c.future;
    }

    Completer c = new Completer();
    resources[path] = c.future;

    loadFile(path).then((bytes) {
      if (bytes == null) {
        c.complete(null);
        return;
      }

      if (_isCompressed(path)) {
        bytes = _decompress(path, bytes);
        _decompressed[path] = true;
      }

      resources[path] = bytes;
      c.complete(bytes);
    });

    return c.future;
  }

  /**
   * Request an image from a file.
   *
   * Requesting the same file multiple times will only load the file once.
   * If the requester wants to process the contents of the file prior to
   * rendering starting, they must pass a [future] to the request. Once the
   * requester recieves the response from the request, it can then process
   * the data and, once finished, complete the [Completer] that the [future]
   * came from. When all requesters of data have completed, the render will
   * continue.
   */
  Future<SpectrumImage> requestImage(String path, [Future future]) {
    if (future != null) {
      futures.add(future);
    }

    if (resources.containsKey(path)) {
      if (resources[path] is Future) {
        return resources[path];
      }

      // If the resource is in raw bytes, then it was probably loaded from an
      // archive and we should decode it now that it's been requested.
      if (resources[path] is List<int>) {
        _decodeImage(path, resources[path]);
      }

      Completer<List<int>> c = new Completer<List<int>>();
      c.complete(resources[path]);

      return c.future;
    }

    LogDebug('LOADING $path');
    Completer<SpectrumImage> c = new Completer<SpectrumImage>();
    resources[path] = c.future;

    loadFile(path).then((bytes) {
      if (bytes == null) {
        LogInfo('UNABLE TO LOAD $path');
        c.complete(null);
        return;
      }
      c.complete(_decodeImage(path, bytes));
    });

    return c.future;
  }

  SpectrumImage _decodeImage(String path, List<int> bytes) {
    Img.Decoder decoder = Img.findDecoderForData(bytes);

    if (decoder == null) {
      resources[path] = null;
      return null;
    }

    Img.HdrImage hdr = decoder.decodeHdrImage(bytes);
    SpectrumImage res = new SpectrumImage(hdr.width, hdr.height);

    int ri = 0;
    for (int y = 0; y < hdr.height; ++y) {
      for (int x = 0; x < hdr.width; ++x) {
        double r = hdr.getRed(x, y);
        double g = hdr.getGreen(x, y);
        double b = hdr.getBlue(x, y);
        res.data[ri++] = r;
        res.data[ri++] = g;
        res.data[ri++] = b;
      }
    }

    LogDebug('IMAGE LOADED $path');
    resources[path] = res;

    return res;
  }

  /**
   * Wait until all requesters of resources have completed their preprocessing.
   */
  Future waitUntilReady() {
    Completer c = new Completer();
    Future.wait(futures).then((r) {
      futures.clear();
      c.complete();
    });
    return c.future;
  }

  /**
   * Get the contents of an already loaded resource.
   */
  getResource(String path) {
    if (!resources.containsKey(path)) {
      return null;
    }

    var data = resources[path];
    if (_isCompressed(path)) {
      data = _decompress(path, data);
      _decompressed[path] = true;
      resources[path] = data;
    }

    return data;
  }

  /**
   * Add a resource, or replace a resources data.
   */
  void setResource(String path, data) {
    resources[path] = data;
  }

  /**
   * Add a resource so it can be accessed later.
   */
  void writeFile(String path, data) {
    resources[path] = data;
  }

  bool hasTexture(String name) => textures.containsKey(name);

  MIPMap getTexture(String name) {
    if (textures.containsKey(name)) {
      return textures[name];
    }
    return null;
  }

  void addTexture(String name, MIPMap texture) {
    textures[name] = texture;
  }

  List<int> _decompress(String file, List<int> bytes) {
    if (file.endsWith('.gz')) {
      bytes = new GZipDecoder().decodeBytes(bytes);
    } else if (file.endsWith('.z')) {
      bytes = new ZLibDecoder().decodeBytes(bytes);
    } else if (file.endsWith('.bz2')) {
      bytes = new BZip2Decoder().decodeBytes(bytes);
    }
    return bytes;
  }

  bool _isCompressed(String file) {
    if (_decompressed.containsKey(file)) {
      if (_decompressed[file]) {
        return false;
      }
    }

    return file.endsWith('.gz') ||
           file.endsWith('.z') ||
           file.endsWith('.bz2');
  }

  List<Future> futures = [];
  Map<String, dynamic> resources = {};
  Map<String, bool> _decompressed = {};
  Map<String, MIPMap> textures = {};
}
