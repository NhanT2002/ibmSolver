use IO;
import input.inputsConfig;
use Time;
use mesh;
use writeCGNS;
use spatialDiscretization;

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

    var FVM = new spatialDiscretization(mesh, inputs);
    FVM.initializeFlowField();
    FVM.updateGhostCells();
    FVM.compute_convective_fluxes();
    FVM.compute_lambdas();

    var writer = new cgnsFlowWriter_c(inputs.OUTPUT_FILENAME_);
    writer.writeToCGNS(mesh);

    time.stop();
    writeln("Runtime: ", time.elapsed(), " seconds");
}