module mesh {

use HDF5;
use HDF5.C_HDF5;

// === Read a 2D dataset from HDF5 ===
proc read2DDataset(type eltType, file_id: hid_t, name: string) {
    // Open the dataset
    var dset_id = H5Dopen(file_id, name.c_str(), H5P_DEFAULT);
    if dset_id < 0 then halt("Dataset not found: ", name);

    // Get dataspace and dimensions
    var space_id = H5Dget_space(dset_id);
    var ndims = H5Sget_simple_extent_ndims(space_id);
    if ndims != 2 then halt("Dataset ", name, " is not 2D");

    var dims: [0..<2] uint(64);
    H5Sget_simple_extent_dims(space_id, c_ptrTo(dims), nil);
    writeln("Reading ", name, " dims = ", dims);

    // Build 2D domain
    var dom = {1..dims[0], 1..dims[1]};
    var data: [dom] eltType;

    // Read dataset
    readHDF5Dataset(file_id, name, data);

    H5Sclose(space_id);
    H5Dclose(dset_id);
    return data;
}

// === Read mesh coordinates X, Y, Z ===
proc readMesh(filename: string) {
    const dsetX = "/Base/Block_1/GridCoordinates/CoordinateX/ data";
    const dsetY = "/Base/Block_1/GridCoordinates/CoordinateY/ data";
    const dsetZ = "/Base/Block_1/GridCoordinates/CoordinateZ/ data";

    var file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if file_id < 0 then halt("Could not open file: ", filename);

    var X = read2DDataset(real(64), file_id, dsetX);
    var Y = read2DDataset(real(64), file_id, dsetY);
    var Z = read2DDataset(real(64), file_id, dsetZ);

    H5Fclose(file_id);
    return (X, Y, Z);
}

class meshData {
    var domain_ : domain(2) = {1..0, 1..0};
    var X_ : [domain_] real(64);
    var Y_ : [domain_] real(64);
    var Z_ : [domain_] real(64);

    var niNode_ : int;
    var njNode_ : int;
    var niCell_ : int;
    var njCell_ : int;

    var node_dom_ : domain(1) = {1..0};
    var xNodes_ : [node_dom_] real(64);
    var yNodes_ : [node_dom_] real(64);
    var zNodes_ : [node_dom_] real(64);

    var cell_dom_ : domain(1) = {1..0};
    var xCells_ : [cell_dom_] real(64);
    var yCells_ : [cell_dom_] real(64);
    var zCells_ : [cell_dom_] real(64);

    var face_dom_ : domain(1) = {1..0};
    var IfaceAreas_ : [face_dom_] real(64);
    var JfaceAreas_ : [face_dom_] real(64);

    proc init(X : [] real(64), Y : [] real(64), Z : [] real(64)) {

        this.domain_ = {1..X.dim(0).size, 1..X.dim(1).size};
        this.X_ = X;
        this.Y_ = Y;
        this.Z_ = Z;

        // flatten the 2D arrays into 1D arrays
        this.niNode_ = X.dim(0).size;
        this.njNode_ = X.dim(1).size;
        this.niCell_ = niNode_ - 1;
        this.njCell_ = njNode_ - 1;

        writeln("Mesh dimensions: niNode = ", niNode_, ", njNode = ", njNode_);

        this.node_dom_ = {1..(niNode_ * njNode_)};
        this.cell_dom_ = {1..(niCell_ * njCell_)};
    }

    proc computeMetrics() {

        for j in 1..njNode_ {
            for i in 1..niNode_ {
                var idx = iiNode(i, j);
                this.xNodes_[idx] = this.X_[j, i];
                this.yNodes_[idx] = this.Y_[j, i];
                this.zNodes_[idx] = this.Z_[j, i];
            }
        }
        for j in 1..njCell_ {
            for i in 1..niCell_ {
                var idx = iiCell(i, j);
                this.xCells_[idx] = 0.25*(this.X_[j, i] + this.X_[j, i+1] + this.X_[j+1, i] + this.X_[j+1, i+1]);
                this.yCells_[idx] = 0.25*(this.Y_[j, i] + this.Y_[j, i+1] + this.Y_[j+1, i] + this.Y_[j+1, i+1]);
                this.zCells_[idx] = 0.25*(this.Z_[j, i] + this.Z_[j, i+1] + this.Z_[j+1, i] + this.Z_[j+1, i+1]);
            }
        }

        // writeln("xCells_ = ", this.xCells_);
        // writeln("yCells_ = ", this.yCells_);
    }

    proc iiNode(i : int, j : int) : int {
        return i + (j-1) * this.niNode_;
    }

    proc iiCell(i : int, j : int) : int {
        return i + (j-1) * this.niCell_;
    }
}



}
