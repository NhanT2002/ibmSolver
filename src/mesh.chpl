module mesh {

use HDF5;
use HDF5.C_HDF5;
use List;
use Set;
use Sort;
use kExactLeastSquare;

// === Generic 1D reader (for real(64) or int) ===
proc read1DDataset(type eltType, file_id: hid_t, name: string): [] eltType {
    var dset_id = H5Dopen(file_id, name.c_str(), H5P_DEFAULT);
    if dset_id < 0 then halt("Dataset not found: ", name);

    var space_id = H5Dget_space(dset_id);
    var dims: [0..<1] uint(64);
    H5Sget_simple_extent_dims(space_id, c_ptrTo(dims), nil);

    var dom = {0..<dims[0]};
    var data: [dom] eltType;
    readHDF5Dataset(file_id, name, data);

    H5Sclose(space_id);
    H5Dclose(dset_id);
    return data;
}

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

    // Build 2D domain
    var dom = {0..<dims[0], 0..<dims[1]};
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

proc readGeometry(filename: string) {
    const dsetX = "/Base/Geometry/GridCoordinates/CoordinateX/ data";
    const dsetY = "/Base/Geometry/GridCoordinates/CoordinateY/ data";
    const dsetZ = "/Base/Geometry/GridCoordinates/CoordinateZ/ data";

    var file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if file_id < 0 then halt("Could not open file: ", filename);
    var X = read1DDataset(real(64), file_id, dsetX);
    var Y = read1DDataset(real(64), file_id, dsetY);
    var Z = read1DDataset(real(64), file_id, dsetZ);

    H5Fclose(file_id);
    return (X, Y, Z);
}

class meshData {
    var geometry_domain_ : domain(1) = {1..0};
    var X_geo_ : [geometry_domain_] real(64);
    var Y_geo_ : [geometry_domain_] real(64);
    var Z_geo_ : [geometry_domain_] real(64);

    var domain_ : domain(2) = {1..0, 1..0};
    var X_ : [domain_] real(64);
    var Y_ : [domain_] real(64);
    var Z_ : [domain_] real(64);

    var niNode_ : int;
    var njNode_ : int;
    var niCell_ : int;
    var njCell_ : int;
    var nCell_ : int;

    var node_dom_ : domain(1) = {1..0};
    var xNodes_ : [node_dom_] real(64);
    var yNodes_ : [node_dom_] real(64);
    var zNodes_ : [node_dom_] real(64);

    var cell_dom_ : domain(1) = {1..0};
    var xCells_ : [cell_dom_] real(64);
    var yCells_ : [cell_dom_] real(64);
    var zCells_ : [cell_dom_] real(64);
    var phiCells_ : [cell_dom_] real(64);
    var cellTypes_ : [cell_dom_] int; // -1: solid, 1: fluid, 0: ghost
    var cellVolumes_ : [cell_dom_] real(64);

    var face_dom_ : domain(1) = {1..0};
    var IfaceAreas_ : [face_dom_] real(64);
    var IfacesCx_ : [face_dom_] real(64);
    var IfacesCy_ : [face_dom_] real(64);
    var JfaceAreas_ : [face_dom_] real(64);
    var JfacesCx_ : [face_dom_] real(64);
    var JfacesCy_ : [face_dom_] real(64);

    var ghostCellDom_ : domain(1) = {1..0};
    var ghostCellIndices_ : [ghostCellDom_] int;
    var ghostCellIJ_ : [ghostCellDom_] (int, int);
    var ghostCells_nx_bi_ : [ghostCellDom_] real(64);
    var ghostCells_ny_bi_ : [ghostCellDom_] real(64);
    var ghostCells_x_mirror_ : [ghostCellDom_] real(64);
    var ghostCells_y_mirror_ : [ghostCellDom_] real(64);
    var ghostCells_x_bi_ : [ghostCellDom_] real(64);
    var ghostCells_y_bi_ : [ghostCellDom_] real(64);

    var ghostCellsNearestFluidCells_dom_ : domain(2) = {1..0, 1..0};
    var ghostCellsNearestFluidCellsCx_dom_ : domain(2) = {1..0, 1..0};
    var ghostCellsNearestFluidCells_ : [ghostCellsNearestFluidCells_dom_] int;
    var ghostCellsNearestFluidCellsCx_ : [ghostCellsNearestFluidCellsCx_dom_] real(64);
    var ghostCellsNearestFluidCellsCy_ : [ghostCellsNearestFluidCellsCx_dom_] real(64);
    var ghostCellkls_ : [ghostCellDom_] owned kls?;


    proc init(X : [] real(64), Y : [] real(64), Z : [] real(64)) {

        this.domain_ = {0..<X.dim(0).size, 0..<X.dim(1).size};
        this.X_ = X;
        this.Y_ = Y;
        this.Z_ = Z;

        // flatten the 2D arrays into 1D arrays
        this.niNode_ = X.dim(0).size;
        this.njNode_ = X.dim(1).size;
        this.niCell_ = niNode_ - 1;
        this.njCell_ = njNode_ - 1;
        this.nCell_ = niCell_ * njCell_;

        writeln("Mesh dimensions: niNode = ", niNode_, ", njNode = ", njNode_);

        this.node_dom_ = {0..<(niNode_ * njNode_)};
        this.cell_dom_ = {0..<(niCell_ * njCell_)};
        this.face_dom_ = {0..<( (niCell_+1)*njCell_ )};
    }

    proc computeMetrics() {

        forall j in 0..<njNode_ {
            for i in 0..<niNode_ {
                var idx = iiNode(i, j);
                this.xNodes_[idx] = this.X_[j, i];
                this.yNodes_[idx] = this.Y_[j, i];
                this.zNodes_[idx] = this.Z_[j, i];
            }
        }

        // Compute cell centroids
        forall j in 0..<njCell_ {
            for i in 0..<niCell_ {
                var idx = iiCell(i, j);
                this.xCells_[idx] = 0.25*(this.X_[j, i] + this.X_[j, i+1] + this.X_[j+1, i] + this.X_[j+1, i+1]);
                this.yCells_[idx] = 0.25*(this.Y_[j, i] + this.Y_[j, i+1] + this.Y_[j+1, i] + this.Y_[j+1, i+1]);
                this.zCells_[idx] = 0.25*(this.Z_[j, i] + this.Z_[j, i+1] + this.Z_[j+1, i] + this.Z_[j+1, i+1]);
            }
        }

        // Compute face areas and centroids (2D)
        forall j in 0..<njCell_ {
            for i in 0..<(niCell_+1) {
                var faceIdx = i + j*(this.niNode_);

                // I-faces (vertical)
                var x0 = this.X_[j, i];
                var y0 = this.Y_[j, i];
                var z0 = this.Z_[j, i];
                var x1 = this.X_[j+1, i];
                var y1 = this.Y_[j+1, i];
                var z1 = this.Z_[j+1, i];
                var dy = y1 - y0;
                this.IfaceAreas_[faceIdx] = abs(dy);
                this.IfacesCx_[faceIdx] = 0.5 * (x0 + x1);
                this.IfacesCy_[faceIdx] = 0.5 * (y0 + y1);
            }
        }

        forall j in 0..<(njCell_+1) {
            for i in 0..<(niCell_) {
                var faceIdx = i + j * (this.niCell_);
                
                // J-faces (horizontal)
                var x0 = this.X_[j, i];
                var y0 = this.Y_[j, i];
                var z0 = this.Z_[j, i];
                var x1 = this.X_[j, i+1];
                var y1 = this.Y_[j, i+1];
                var z1 = this.Z_[j, i+1];
                var dx = x1 - x0;
                this.JfaceAreas_[faceIdx] = abs(dx);
                this.JfacesCx_[faceIdx] = 0.5 * (x0 + x1);
                this.JfacesCy_[faceIdx] = 0.5 * (y0 + y1);
            }
        }

        // Compute cell volumes
        forall j in 0..<njCell_ {
            for i in 0..<niCell_ {
                var idx = iiCell(i, j);
                var x1 = this.X_[j, i];
                var y1 = this.Y_[j, i];
                var x2 = this.X_[j, i+1];
                var y2 = this.Y_[j, i+1];
                var x3 = this.X_[j+1, i+1];
                var y3 = this.Y_[j+1, i+1];
                var x4 = this.X_[j+1, i];
                var y4 = this.Y_[j+1, i];
                this.cellVolumes_[idx] = 0.5 * ((x1-x3)*(y2-y4) + (x4-x2)*(y1-y3));
            }
        }
    }

    proc iiNode(i : int, j : int) : int {
        return i + j * this.niNode_;
    }

    proc iiCell(i : int, j : int) : int {
        return i + j * this.niCell_;
    }

    proc levelSet(X_geo : [] real(64), Y_geo : [] real(64)) {
        this.geometry_domain_ = {0..<X_geo.size};
        this.X_geo_ = X_geo;
        this.Y_geo_ = Y_geo;
        // build list of segments (A->B) with wrap around
        record Seg {
            var Ax: real(64);
            var Ay: real(64);
            var Bx: real(64);
            var By: real(64);
        }
        var segs : [{0..<X_geo.size}] Seg;
        for k in 0..<X_geo.size {
            var next = (k + 1) % X_geo.size;
            segs[k] = new Seg(Ax = X_geo[k], Ay = Y_geo[k],
                             Bx = X_geo[next], By = Y_geo[next]);
        }

        // helper function to compute distance from point to segment
        proc pointToSegmentDistance(px: real(64), py: real(64), const ref seg: Seg) : real(64) {
            var Ax = seg.Ax;
            var Ay = seg.Ay;
            var Bx = seg.Bx;
            var By = seg.By;

            var ABx = Bx - Ax;
            var ABy = By - Ay;
            var APx = px - Ax;
            var APy = py - Ay;
            var AB_squared = ABx*ABx + ABy*ABy;
            var t = 0.0;
            if AB_squared > 0.0 {
                t = (APx*ABx + APy*ABy) / AB_squared;
            }

            var cx = 0.0;
            var cy = 0.0;
            if t <= 0.0 {
                cx = Ax;
                cy = Ay;
            } else if t >= 1.0 {
                cx = Bx;
                cy = By;
            } else {
                cx = Ax + t * ABx;
                cy = Ay + t * ABy;
            }
            var dx = px - cx;
            var dy = py - cy;
            return sqrt(dx*dx + dy*dy);
        }

        // helper: test horizontal ray intersection with segment
        proc rayIntersectsSegment(px: real(64), py: real(64), const ref seg: Seg) : bool {
            var Ay = seg.Ay;
            var By = seg.By;
            var Ax = seg.Ax;
            var Bx = seg.Bx;
            if (Ay == By) then return false; // horizontal segment
            var cond = ((Ay > py) != (By > py));
            if !cond then return false;
            var xint = Ax + (py - Ay) * (Bx - Ax) / (By - Ay);
            return (xint > px);
        }

        // loop over cells to compute level set
        forall idx in 0..<nCell_ {
            var px = xCells_[idx];
            var py = yCells_[idx];

            var minDist = 1.0e20;
            for s in segs {
                var dist = pointToSegmentDistance(px, py, s);
                if dist < minDist {
                    minDist = dist;
                }
            }

            var intersections : int = 0;
            for s in segs {
                if rayIntersectsSegment(px, py, s) {
                    intersections += 1;
                }
            }

            if (intersections % 2) == 1 {
                // inside geometry
                this.phiCells_[idx] = -minDist;
                this.cellTypes_[idx] = -1;
            } else {
                // outside geometry
                this.phiCells_[idx] = minDist;
                this.cellTypes_[idx] = 1;
            }

        }
        

        // Identify ghost cells
        var ghostCellList = new list((int, int));
        for j in 0..<njCell_ {
            for i in 0..<niCell_ {
                var idx = iiCell(i, j);
                if cellTypes_[idx] == -1 {
                    // solid cell, check neighbors
                    var neighbors = [idx - 1, idx + 1, idx - niCell_, idx + niCell_];
                    for nidx in neighbors {
                        if this.cellTypes_[nidx] == 1 {
                            // neighbor is fluid, mark as ghost
                            this.cellTypes_[idx] = 0;
                            ghostCellList.pushBack((i, j));
                            break;
                        }
                    }
                }
            }
        }

        this.ghostCellDom_ = {0..<ghostCellList.size};
        forall i in this.ghostCellDom_ {
            this.ghostCellIJ_[i] = ghostCellList[i];
            this.ghostCellIndices_[i] = iiCell(ghostCellIJ_[i][0], ghostCellIJ_[i][1]);
        }

        this.ghostCellsNearestFluidCells_dom_ = {0..<ghostCellDom_.size, 0..2};
        this.ghostCellsNearestFluidCellsCx_dom_ = {0..<ghostCellDom_.size, 0..3};
    }

    proc computeIBnormals() {
        record Seg {
            var Ax: real(64);
            var Ay: real(64);
            var Bx: real(64);
            var By: real(64);
        }
        var segs : [{0..<this.X_geo_.size}] Seg;
        for k in 0..<this.X_geo_.size {
            var next = (k + 1) % this.X_geo_.size;
            segs[k] = new Seg(Ax = this.X_geo_[k], Ay = this.Y_geo_[k],
                             Bx = this.X_geo_[next], By = this.Y_geo_[next]);
        }

        // helper function to compute projection of point onto segment
        proc project_to_seg(px: real(64), py: real(64), const ref seg: Seg, ref cxp: real(64), ref cyp: real(64)) {
            const ax = seg.Ax;
            const ay = seg.Ay;
            const bx = seg.Bx;
            const By = seg.By;
            const abx = bx - ax;
            const aby = By - ay;
            const apx = px - ax;
            const apy = py - ay;
            const ab2 = abx*abx + aby*aby;
            var t = 0.0;
            if ab2 > 0.0 {
                t = (apx*abx + apy*aby) / ab2;
            }
            if t <= 0.0 {
                cxp = ax;
                cyp = ay;
            } else if t > 1.0 {
                cxp = bx;
                cyp = By;
            } else {
                cxp = ax + t * abx;
                cyp = ay + t * aby;
            }
            const dx = px - cxp;
            const dy = py - cyp;
            const dist2 = dx*dx + dy*dy;
            return sqrt(dist2);
        }

        // loop over ghost cells to compute normals
        forall i in 0..<ghostCellDom_.size {
            const ghostCell = ghostCellIndices_[i];
            var ox = xCells_[ghostCell];
            var oy = yCells_[ghostCell];
            var min_d = 1.0e20;
            var best_x = ox;
            var best_y = oy;
            for s in segs {
                var cxp = 0.0;
                var cyp = 0.0;
                var d = project_to_seg(ox, oy, s, cxp, cyp);
                if d < min_d {
                    min_d = d;
                    best_x = cxp;
                    best_y = cyp;
                }
            }
            const vx = best_x - ox;
            const vy = best_y - oy;
            const norm = sqrt(vx*vx + vy*vy);
            const nx = vx / norm;
            const ny = vy / norm;
            this.ghostCells_nx_bi_[i] = nx;
            this.ghostCells_ny_bi_[i] = ny;
            this.ghostCells_x_mirror_[i] = best_x + vx;
            this.ghostCells_y_mirror_[i] = best_y + vy;
            this.ghostCells_x_bi_[i] = best_x;
            this.ghostCells_y_bi_[i] = best_y;

            // Find nearest 3 fluid cell
            var dists = new list((real(64), int));
            var potentialCells = new set(int);
            // neighboring fluid cells
            var neighbors = [ghostCell-1, ghostCell+1, ghostCell-niCell_, ghostCell+niCell_];
            for nidx in neighbors {
                potentialCells.add(nidx);
                // also add their neighbors
                var neighborNeighbors = [nidx-1, nidx+1, nidx-niCell_, nidx+niCell_];
                for nnidx in neighborNeighbors {
                    potentialCells.add(nnidx);
                    // also add their neighbors
                    var nnNeighbors = [nnidx-1, nnidx+1, nnidx-niCell_, nnidx+niCell_];
                    for nnnidx in nnNeighbors {
                        potentialCells.add(nnnidx);
                    }
                }
            }
            for cell in potentialCells {
                if cellTypes_[cell] == 1 {
                    var dx = xCells_[cell] - this.ghostCells_x_mirror_[i];
                    var dy = yCells_[cell] - this.ghostCells_y_mirror_[i];
                        var dist = sqrt(dx*dx + dy*dy);
                        dists.pushBack((dist, cell));
                    }
                }
            // Sort distances and get indices of nearest 3
            sort(dists);
            for k in 0..<3 {
                this.ghostCellsNearestFluidCells_[i, k] = dists[k][1];
                this.ghostCellsNearestFluidCellsCx_[i, k] = xCells_[dists[k][1]];
                this.ghostCellsNearestFluidCellsCy_[i, k] = yCells_[dists[k][1]];
            }
            this.ghostCellsNearestFluidCellsCx_[i, 3] = this.ghostCells_x_bi_[i];
            this.ghostCellsNearestFluidCellsCy_[i, 3] = this.ghostCells_y_bi_[i];
            this.ghostCellkls_[i] = new owned kls(this.ghostCellsNearestFluidCellsCx_[i, 0..],
                                                  this.ghostCellsNearestFluidCellsCy_[i, 0..],
                                                  this.ghostCells_x_mirror_[i],
                                                  this.ghostCells_y_mirror_[i]);
            this.ghostCellkls_[i]!.computeCoefficients();
        }
        
    }

}

} // module mesh
