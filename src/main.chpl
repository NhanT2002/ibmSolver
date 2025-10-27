use IO;
import input.inputsConfig;
use Time;
use mesh;
use cgns;

proc main() {
    var time: stopwatch;
    time.start();

    var inputs = new inputsConfig();

    var (X, Y, Z) = readMesh(inputs.MESH_FILENAME_);

    var mesh = new meshData(X, Y, Z);
    mesh.computeMetrics();

    time.stop();
    writeln("Runtime: ", time.elapsed(), " seconds");
}