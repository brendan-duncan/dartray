part of dartray_web;

class RenderManager extends RenderManagerInterface {
  Future<List<int>> requrestFile(String path) {
    print('LOAD $path');
    Completer<List<int>> c = new Completer<List<int>>();
    c.complete(null);
    return c.future;
  }
}
