use IO;
import input.inputsConfig;
use Time;
use mesh;
use spatialDiscretization;
use temporalDiscretization;
use fullPotentialSpatialDiscretization;
use fullPotentialTemporalDiscretization;
use List;

proc main() {
    var time: stopwatch;
    time.start();

    var inputs = new inputsConfig();
    inputs.initializeFlowField();

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
        var time_list = new list(real);
        var iterations_list = new list(int);
        var res0_list = new list(real);
        var res1_list = new list(real);
        var cls_list = new list(real);
        var cds_list = new list(real);
        var cms_list = new list(real);
        time_list.pushBack(0.0);
        iterations_list.pushBack(0);
        res0_list.pushBack(0.0);
        res1_list.pushBack(0.0);
        cls_list.pushBack(0.0);
        cds_list.pushBack(0.0);
        cms_list.pushBack(0.0);
        FVM.writeSolution2CGNS(time_list, iterations_list, res0_list, res1_list, cls_list, cds_list, cms_list);

        // var solver = new shared fullPotentialTemporalDiscretization(FVM, inputs);
        // solver.solve();
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