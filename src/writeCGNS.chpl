module writeCGNS 
{
use cgns;
use ElementTopology;
use CTypes;
use mesh;
use FileSystem;
use List;
use Map;
use arrayOperations;
use CGNSextern;

class MyCGNSfile_c : CGNSfile_c {
    proc init(filename: string, mode: CGNSOpenMode_t, basename: string) {
        super.init(filename, mode, basename);
    }
}

class cgnsFlowWriter_c {
    var cgnsFileName_ : string;
    var cgnsFile_    : owned MyCGNSfile_c;
    var glbBaseId : c_int;
    var zoneId : c_int;
    var solIDcc : c_int;

    proc init(inputFileName : string) {
        var baseName : string = "Base";

        cgnsFileName_ = inputFileName;

        var name = new list(string);
        for st in inputFileName.split(".") {
        name.pushBack(st);
        }

        var i: int = 0;
        var fileExist = try! exists(cgnsFileName_);
        while fileExist {
        i += 1;
        cgnsFileName_ = name[0] + "_" + i:string + "." + name[1];
        fileExist = try! exists(cgnsFileName_);
        }

        cgnsFile_ = new owned MyCGNSfile_c(cgnsFileName_, CGNSOpenMode_t.WRITE, baseName);

        const cellDim : c_int     = 2;
        const physDim : c_int     = 3;

        glbBaseId = cgnsFile_.createBase(baseName, cellDim, physDim);
    }

    proc writeToCGNS(meshData_ : meshData) {
        const zoneName : string = "Zone";

        var size : [0..5] cgsize_t;
        size[0] = meshData_.niNode_;
        size[1] = meshData_.njNode_;
        size[2] = meshData_.niNode_ - 1;
        size[3] = meshData_.njNode_ - 1;
        size[4] = 0;
        size[5] = 0;

        zoneId = cgnsFile_.createZone(glbBaseId, zoneName, size, Structured);

        cgnsFile_.writeGridCoordinates(glbBaseId, zoneId, meshData_.xNodes_, meshData_.yNodes_, meshData_.zNodes_);

        solIDcc = cgnsFile_.addCellCenteredSolution(glbBaseId, zoneId, "FLOW_SOLUTION_CC");

        var idCells : [meshData_.cell_dom_] int;
        forall i in meshData_.cell_dom_ {
            idCells[i] = i;
        }

        var fields = new map(string, [meshData_.cell_dom_] real(64));
        fields["idCells"] = idCells;
        fields["phiCells"] = meshData_.phiCells_;
        fields["cellTypes"] = meshData_.cellTypes_;
        fields["cellVolumes"] = meshData_.cellVolumes_;
        for name in fields.keys() {
            var values = try! fields[name];
            cgnsFile_.addFieldSolution(glbBaseId, zoneId, solIDcc, name, values);
        }


        cgnsFile_.close();
        writeln("CGNS file written to ", cgnsFileName_);
    }

    proc writeToCGNS(meshData_ : meshData, dom: domain(1), fieldMap: map(string, [dom] real(64))) {
        const zoneName : string = "Zone";

        var size : [0..5] cgsize_t;
        size[0] = meshData_.niNode_;
        size[1] = meshData_.njNode_;
        size[2] = meshData_.niNode_ - 1;
        size[3] = meshData_.njNode_ - 1;
        size[4] = 0;
        size[5] = 0;

        zoneId = cgnsFile_.createZone(glbBaseId, zoneName, size, Structured);

        cgnsFile_.writeGridCoordinates(glbBaseId, zoneId, meshData_.xNodes_, meshData_.yNodes_, meshData_.zNodes_);

        solIDcc = cgnsFile_.addCellCenteredSolution(glbBaseId, zoneId, "FLOW_SOLUTION_CC");

        var idCells : [meshData_.cell_dom_] int;
        forall i in meshData_.cell_dom_ {
            idCells[i] = i;
        }

        var fields = new map(string, [meshData_.cell_dom_] real(64));
        fields["idCells"] = idCells;
        fields["phiCells"] = meshData_.phiCells_;
        fields["cellTypes"] = meshData_.cellTypes_;
        fields["cellVolumes"] = meshData_.cellVolumes_;
        for name in fields.keys() {
            var values = try! fields[name];
            cgnsFile_.addFieldSolution(glbBaseId, zoneId, solIDcc, name, values);
        }

        for name in fieldMap.keys() {
            var values = try! fieldMap[name];
            cgnsFile_.addFieldSolution(glbBaseId, zoneId, solIDcc, name, values);
        }

        writeln("Cell-centered CGNS file written to ", cgnsFileName_);
    }

    proc writeWallToCGNS(meshData_ : meshData, dom: domain(1), fieldMap: map(string, [dom] real(64))) {

        const wallBaseName = "WallBase";
        const cellDim: c_int = 1;      // 1D elements (lines)
        const physDim: c_int = 2;      // 2D space (x, y)
        const wallBaseId = cgnsFile_.createBase(wallBaseName, cellDim, physDim);

        const zoneName : string = "ZoneWall";
        var size : [0..2] cgsize_t;
        size[0] = meshData_.ghostCellDom_.size;
        size[1] = meshData_.ghostCellDom_.size-1;
        size[2] = 0;

        zoneId = cgnsFile_.createZone(wallBaseId, zoneName, size, Structured);

        cgnsFile_.writeGridCoordinates(wallBaseId, zoneId, meshData_.ghostCells_x_bi_, meshData_.ghostCells_y_bi_, meshData_.ghostCells_z_bi_);

        solIDcc = cgnsFile_.addNodeCenteredSolution(wallBaseId, zoneId, "WALL_FLOW_SOLUTION_NC");

        for name in fieldMap.keys() {
            var values = try! fieldMap[name];
            cgnsFile_.addFieldSolution(wallBaseId, zoneId, solIDcc, name, values);
        }

        writeln("CGNS file written to ", cgnsFileName_);


        cgnsFile_.close();
    }

    proc writeConvergenceHistory(time: list(real(64)),
                                 iterations: list(int),
                                 res0: list(real(64)),
                                 res1: list(real(64)),
                                 res2: list(real(64)),
                                 res3: list(real(64)),
                                 cls: list(real(64)),
                                 cds: list(real(64)),
                                 cms: list(real(64))) {
        const nMetrics = 9;
        const maxIter = iterations.size;

        var iterationsArray: [0..<maxIter] int;
        forall i in 0..<maxIter {
            iterationsArray[i] = iterations[i];
        }

        var convData: [0..<nMetrics, 0..<maxIter] real(64);
        convData[0, ..] = time;
        convData[1, ..] = res0;
        convData[2, ..] = res1;
        convData[3, ..] = res2;
        convData[4, ..] = res3;
        convData[5, ..] = cls;
        convData[6, ..] = cds;
        convData[7, ..] = cms;

        var names = ["Time", "Res0", "Res1", "Res2", "Res3", "Cl", "Cd", "Cm"];

        cgnsFile_.writeGlobalConvergenceHistory("ConvergenceHistory", iterationsArray, names, convData);

        writeln("Convergence history written to ", cgnsFileName_);

    }
}


} // module writeCGNS