part of dartray_io;

class RenderManager extends RenderManagerInterface {
  Future<OutputImage> renderFile(String path) {
    File fp = new File(path);
    addIncludePath(fp.parent.path);
    String scene = fp.readAsStringSync();
    return render(scene);
  }

  Future<List<int>> requrestFile(String file) {
    String path = _findFile(file);
    Completer<List<int>> c = new Completer<List<int>>();

    if (path == null) {
      LogInfo('File not found: $file');
      c.complete(null);
      return c.future;
    }

    LogInfo('LOAD $path');
    Future<List<int>> f = new File(path).readAsBytes();
    f.then((bytes) {
      resources[file] = bytes;
    });

    return f;
  }

  List<int> getFile(String path) {
    if (!resources.containsKey(path)) {
      return null;
    }
    return resources[path];
  }

  String _findFile(String file) {
    if (FileSystemEntity.isFileSync(file)) {
      return file;
    }
    for (String p in includePaths) {
      if (FileSystemEntity.isFileSync(p + '/' + file)) {
        return p + '/' + file;
      }
    }
    return null;
  }
}
