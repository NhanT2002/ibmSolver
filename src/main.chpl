use IO;
import input.inputsConfig;
use Time;
use mesh;
use spatialDiscretization;
use temporalDiscretization;

proc main() {
    var time: stopwatch;
    time.start();

    var inputs = new inputsConfig();

    var (X, Y, Z) = readMesh(inputs.MESH_FILENAME_);
    var (X_geo, Y_geo, Z_geo) = readGeometry(inputs.GEOMETRY_FILENAME_);

    var mesh = new shared meshData(X, Y, Z);
    mesh.computeMetrics();
    mesh.levelSet(X_geo, Y_geo);
    mesh.computeIBnormals();

    var FVM = new shared spatialDiscretization(mesh, inputs);
    FVM.initializeFlowField();
    FVM.run_odd();

    var solver = new shared temporalDiscretization(FVM, inputs);
    solver.solve();

    time.stop();
    writeln("Runtime: ", time.elapsed(), " seconds");
}