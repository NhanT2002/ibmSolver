use IO;
import input.inputsConfig;
use Time;
use mesh;
use spatialDiscretization;
use temporalDiscretization;
use fullPotentialSpatialDiscretization;
use fullPotentialTemporalDiscretization;

proc main() {
    var time: stopwatch;
    time.start();

    var inputs = new inputsConfig();

    var (X, Y, Z) = readMesh(inputs.MESH_FILENAME_);
    var (X_geo, Y_geo, Z_geo) = readGeometry(inputs.GEOMETRY_FILENAME_);

    var mesh = new shared meshData(inputs, X, Y, Z);
    mesh.computeMetrics();
    mesh.levelSet(X_geo, Y_geo);
    mesh.computeIBnormals();

    if inputs.FLOW_ == "fullPotential" {
        var FVM = new shared fullPotentialSpatialDiscretization(mesh, inputs);
        FVM.initializeFlowField();
        FVM.run();

        var solver = new shared fullPotentialTemporalDiscretization(FVM, inputs);
        solver.solve();
    } 
    else {
        var FVM = new shared spatialDiscretization(mesh, inputs);
        FVM.initializeFlowField();
        FVM.run_odd();

        var solver = new shared temporalDiscretization(FVM, inputs);
        solver.solve();
    }

    time.stop();
    writeln("Runtime: ", time.elapsed(), " seconds");
}