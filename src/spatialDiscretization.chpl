module spatialDiscretization 
{
use mesh;
use Time;
use Random;
use List;
import input.inputsConfig;

class spatialDiscretization {
    var time_: stopwatch;
    var mesh_ : shared meshData;
    var inputs_ : inputsConfig;

    var cell_dom_with_ghosts_ : domain(1) = {1..0};
    var niCellWithGhosts_ : int;
    var njCellWithGhosts_ : int;
    var xCellsWithGhosts_ : [cell_dom_with_ghosts_] real(64);
    var yCellsWithGhosts_ : [cell_dom_with_ghosts_] real(64);
    var cellTypesWithGhosts_ : [cell_dom_with_ghosts_] int;

    var W0_ : [cell_dom_with_ghosts_] real(64);
    var W1_ : [cell_dom_with_ghosts_] real(64);
    var W2_ : [cell_dom_with_ghosts_] real(64);
    var W3_ : [cell_dom_with_ghosts_] real(64);

    var rhorho_ : [cell_dom_with_ghosts_] real(64);
    var uu_ : [cell_dom_with_ghosts_] real(64);
    var vv_ : [cell_dom_with_ghosts_] real(64);
    var EE_ : [cell_dom_with_ghosts_] real(64);
    var pp_ : [cell_dom_with_ghosts_] real(64);

    var Rc0_ : [cell_dom_with_ghosts_] real(64);
    var Rc1_ : [cell_dom_with_ghosts_] real(64);
    var Rc2_ : [cell_dom_with_ghosts_] real(64);
    var Rc3_ : [cell_dom_with_ghosts_] real(64);

    var LambdaI_ : [cell_dom_with_ghosts_] real(64);
    var LambdaJ_ : [cell_dom_with_ghosts_] real(64);

    var Rd0_ : [cell_dom_with_ghosts_] real(64);
    var Rd1_ : [cell_dom_with_ghosts_] real(64);
    var Rd2_ : [cell_dom_with_ghosts_] real(64);
    var Rd3_ : [cell_dom_with_ghosts_] real(64);

    var ghostCellIndices_dom_ : domain(1) = {1..0};
    var ghostCellWallIndicesWithGhost_ : [ghostCellIndices_dom_] int;
    var ghostCellsNearestFluidCellsWithGhost_ : [ghostCellIndices_dom_] (int, int, int);

    var ghostCells_rho_bi_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_u_bi_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_v_bi_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_E_bi_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_p_bi_ : [ghostCellIndices_dom_] real(64);

    var ghostCells_rho_mirror_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_u_mirror_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_v_mirror_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_E_mirror_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_p_mirror_ : [ghostCellIndices_dom_] real(64);

    var face_dom_ : domain(1) = {1..0};
    var F0I_ : [face_dom_] real(64);
    var F1I_ : [face_dom_] real(64);
    var F2I_ : [face_dom_] real(64);
    var F3I_ : [face_dom_] real(64);

    var F0J_ : [face_dom_] real(64);
    var F1J_ : [face_dom_] real(64);
    var F2J_ : [face_dom_] real(64);
    var F3J_ : [face_dom_] real(64);

    var D0I_ : [face_dom_] real(64);
    var D1I_ : [face_dom_] real(64);
    var D2I_ : [face_dom_] real(64);
    var D3I_ : [face_dom_] real(64);

    var D0J_ : [face_dom_] real(64);
    var D1J_ : [face_dom_] real(64);
    var D2J_ : [face_dom_] real(64);
    var D3J_ : [face_dom_] real(64);

    var outflowFacesI_dom_ : domain(1) = {1..0};
    var inflowFacesI_dom_ : domain(1) = {1..0};
    var outflowFacesJ_dom_ : domain(1) = {1..0};
    var inflowFacesJ_dom_ : domain(1) = {1..0};

    var outflowFacesI_ : [outflowFacesI_dom_] int;
    var outflowFacesJ_ : [outflowFacesJ_dom_] int;
    var inflowFacesI_ : [inflowFacesI_dom_] int;
    var inflowFacesJ_ : [inflowFacesJ_dom_] int;

    proc init(Mesh: shared meshData, ref inputs: inputsConfig) {
        this.time_ = new stopwatch();
        this.mesh_ = Mesh;
        this.inputs_ = inputs;

        // Define domains
        const niCell_ = mesh_.niCell_;
        const njCell_ = mesh_.njCell_;
        this.cell_dom_with_ghosts_ = {0..<( (niCell_ + 4) * (njCell_ + 4) )};
        this.niCellWithGhosts_ = niCell_ + 4;
        this.njCellWithGhosts_ = njCell_ + 4;
        this.ghostCellIndices_dom_ = this.mesh_.ghostCellIndices_.domain;
        this.face_dom_ = {0..<( (niCell_+1)*njCell_ )};

        writeln("spatialDiscretization initialized.");
        writeln("  niCell_ = ", niCell_);
        writeln("  njCell_ = ", njCell_);
        writeln("  niCellWithGhosts_ = ", this.niCellWithGhosts_);
        writeln("  njCellWithGhosts_ = ", this.njCellWithGhosts_);
    }

    proc initializeFlowField() {
        this.W0_ = this.inputs_.RHO_INF_;
        this.W1_ = this.inputs_.RHO_INF_ * this.inputs_.U_INF_;
        this.W2_ = this.inputs_.RHO_INF_ * this.inputs_.V_INF_;
        this.W3_ = this.inputs_.RHO_INF_ * this.inputs_.E_INF_;

        this.rhorho_ = this.inputs_.RHO_INF_;
        this.uu_ = this.inputs_.U_INF_;
        this.vv_ = this.inputs_.V_INF_;
        this.EE_ = this.inputs_.E_INF_;
        this.pp_ = this.inputs_.P_INF_;

        this.ghostCells_rho_bi_ = this.inputs_.RHO_INF_;
        this.ghostCells_u_bi_ = this.inputs_.U_INF_;
        this.ghostCells_v_bi_ = this.inputs_.V_INF_;
        this.ghostCells_E_bi_ = this.inputs_.E_INF_;
        this.ghostCells_p_bi_ = this.inputs_.P_INF_;

        // Map to mesh cells with ghosts
        this.xCellsWithGhosts_ = 1e100;
        this.yCellsWithGhosts_ = 1e100;
        forall j in 0..<this.mesh_.njCell_ {
            for i in 0..<this.mesh_.niCell_ {
                const ii = i + 2;
                const jj = j + 2;
                const idxWithGhost = ii + jj * this.niCellWithGhosts_;
                const idx = this.mesh_.iiCell(i, j);
                this.xCellsWithGhosts_[idxWithGhost] = this.mesh_.xCells_[idx];
                this.yCellsWithGhosts_[idxWithGhost] = this.mesh_.yCells_[idx];
            }
        }
        forall i in 2..<this.niCellWithGhosts_-2 {
            this.xCellsWithGhosts_[i+this.niCellWithGhosts_] = this.xCellsWithGhosts_[i + 2 * this.niCellWithGhosts_];
            this.yCellsWithGhosts_[i+this.niCellWithGhosts_] = this.xCellsWithGhosts_[i + 2 * this.niCellWithGhosts_];
        }
        forall i in 2..<this.niCellWithGhosts_-2 {
            this.xCellsWithGhosts_[i] = this.xCellsWithGhosts_[i + 2 * this.niCellWithGhosts_];
            this.yCellsWithGhosts_[i] = this.xCellsWithGhosts_[i + 2 * this.niCellWithGhosts_];
        }
        forall i in 2..<this.niCellWithGhosts_-2 {
            this.xCellsWithGhosts_[i + (this.njCellWithGhosts_-1) * this.niCellWithGhosts_] = this.xCellsWithGhosts_[i + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_];
            this.yCellsWithGhosts_[i + (this.njCellWithGhosts_-1) * this.niCellWithGhosts_] = this.yCellsWithGhosts_[i + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_];
        }
        forall i in 2..<this.niCellWithGhosts_-2 {
            this.xCellsWithGhosts_[i + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] = this.xCellsWithGhosts_[i + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_];
            this.yCellsWithGhosts_[i + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] = this.yCellsWithGhosts_[i + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_];
        }
        // // test 
        // var rands: [0..<this.niCellWithGhosts_] int;
        // forall i in 0..<this.niCellWithGhosts_ {
        //     rands[i] = i;
        // }
        // writeln("rands before = ", rands);
        // var samp = try! sample(rands, 2);
        // writeln("sample = ", samp);
        // const i = samp[0];
        // const j = samp[1];
        // const idxWithGhost = i + j * this.niCellWithGhosts_;
        // const idx = this.mesh_.iiCell(i-2, j-2);
        // writeln("Test mapping: i=", i, " j=", j, " idxWithGhost=", idxWithGhost, " idx=", idx);
        // writeln(" xCellsWithGhosts_=", this.xCellsWithGhosts_[idxWithGhost], " mesh xCells_=", this.mesh_.xCells_[idx]);
        // writeln(" yCellsWithGhosts_=", this.yCellsWithGhosts_[idxWithGhost], " mesh yCells_=", this.mesh_.yCells_[idx]);

        forall (id,(i, j)) in zip(this.mesh_.ghostCellIJ_.domain, this.mesh_.ghostCellIJ_) {
            const idxWithGhost = i+2 + (j+2) * this.niCellWithGhosts_;
            this.ghostCellWallIndicesWithGhost_[id] = idxWithGhost;
        }
        forall i in this.ghostCellIndices_dom_ {
            for j in this.mesh_.ghostCellsNearestFluidCells_dom_.dim(1) {
                const (fi, fj) = this.mesh_.ghostCellsNearestFluidCellsIJ_[i, j];
                const fi_withGhost = fi + 2;
                const fj_withGhost = fj + 2;
                const idxWithGhost = fi_withGhost + fj_withGhost * this.niCellWithGhosts_;
                this.ghostCellsNearestFluidCellsWithGhost_[i][j] = idxWithGhost;
            }
        }
        this.cellTypesWithGhosts_ = 1;
        forall j in 0..<this.mesh_.njCell_ {
            for i in 0..<this.mesh_.niCell_ {
                const idxWithGhost = meshIndex2FVMindex(i, j);
                const idx = this.mesh_.iiCell(i, j);
                this.cellTypesWithGhosts_[idxWithGhost] = this.mesh_.cellTypes_[idx];
            }
        }

        // find outflow and inflow faces
        var outflowFacesI = new list(int);
        var outflowFacesJ = new list(int);
        var inflowFacesI = new list(int);
        var inflowFacesJ = new list(int);

        if this.inputs_.ALPHA_ == 0.0 {
            for j in 0..<this.mesh_.njCell_ {
                inflowFacesI.pushBack(0 + j * (this.mesh_.niCell_ + 1));
                outflowFacesI.pushBack((this.mesh_.niCell_) + j * (this.mesh_.niCell_ + 1));
            }
            for i in 0..<this.mesh_.niCell_ {
                inflowFacesJ.pushBack(i + 0 * this.mesh_.niCell_);
                inflowFacesJ.pushBack(i + this.mesh_.niCell_ * this.mesh_.njCell_);
            }
        }
        else if this.inputs_.ALPHA_ > 0 {
            for j in 0..<this.mesh_.njCell_ {
                inflowFacesI.pushBack(0 + j * (this.mesh_.niCell_ + 1));
                outflowFacesI.pushBack((this.mesh_.niCell_) + j * (this.mesh_.niCell_ + 1));
            }
            for i in 0..<this.mesh_.niCell_ {
                inflowFacesJ.pushBack(i + 0 * this.mesh_.niCell_);
                outflowFacesJ.pushBack(i + this.mesh_.niCell_ * this.mesh_.njCell_);
            }
        }
        else if this.inputs_.ALPHA_ < 0 {
            for j in 0..<this.mesh_.njCell_ {
                inflowFacesI.pushBack(0 + j * (this.mesh_.niCell_ + 1));
                outflowFacesI.pushBack((this.mesh_.niCell_) + j * (this.mesh_.niCell_ + 1));
            }
            for i in 0..<this.mesh_.niCell_ {
                outflowFacesJ.pushBack(i + 0 * this.mesh_.niCell_);
                inflowFacesJ.pushBack(i + this.mesh_.niCell_ * this.mesh_.njCell_);
            }
        }

        this.outflowFacesI_dom_ = {0..<outflowFacesI.size};
        this.outflowFacesJ_dom_ = {0..<outflowFacesJ.size};
        this.inflowFacesI_dom_ = {0..<inflowFacesI.size};
        this.inflowFacesJ_dom_ = {0..<inflowFacesJ.size};

        forall i in this.outflowFacesI_dom_ {
            this.outflowFacesI_[i] = outflowFacesI[i];
        }
        forall i in this.outflowFacesJ_dom_ {
            this.outflowFacesJ_[i] = outflowFacesJ[i];
        }
        forall i in this.inflowFacesI_dom_ {
            this.inflowFacesI_[i] = inflowFacesI[i];
        }
        forall i in this.inflowFacesJ_dom_ {
            this.inflowFacesJ_[i] = inflowFacesJ[i];
        }
    }

    proc updateGhostCells() {
        // Wall
        forall i in this.ghostCellIndices_dom_ {
            const ghostCell = this.ghostCellWallIndicesWithGhost_[i];
            var kls = this.mesh_.ghostCellkls_[i]!.borrow();
            for (j, nearestCell) in zip(0..2, this.ghostCellsNearestFluidCellsWithGhost_[i]) { // MAYBE GENERALIZE LATER
                kls.rho_fieldValues_[j] = this.rhorho_[nearestCell];
                kls.u_fieldValues_[j] = this.uu_[nearestCell];
                kls.v_fieldValues_[j] = this.vv_[nearestCell];
                kls.E_fieldValues_[j] = this.EE_[nearestCell];
                kls.p_fieldValues_[j] = this.pp_[nearestCell];
            }
            // last value is at boundary interface
            kls.rho_fieldValues_[3] = this.ghostCells_rho_bi_[i];
            kls.u_fieldValues_[3] = this.ghostCells_u_bi_[i];
            kls.v_fieldValues_[3] = this.ghostCells_v_bi_[i];
            kls.E_fieldValues_[3] = this.ghostCells_E_bi_[i];
            kls.p_fieldValues_[3] = this.ghostCells_p_bi_[i];

            const rho_mirror = kls.interpolate(kls.rho_fieldValues_);
            const u_mirror = kls.interpolate(kls.u_fieldValues_);
            const v_mirror = kls.interpolate(kls.v_fieldValues_);
            const E_mirror = kls.interpolate(kls.E_fieldValues_);
            const p_mirror = kls.interpolate(kls.p_fieldValues_);

            this.uu_[ghostCell] = u_mirror - 2.0 * (u_mirror * this.mesh_.ghostCells_nx_bi_[i] + v_mirror * this.mesh_.ghostCells_ny_bi_[i]) * this.mesh_.ghostCells_nx_bi_[i];
            this.vv_[ghostCell] = v_mirror - 2.0 * (u_mirror * this.mesh_.ghostCells_nx_bi_[i] + v_mirror * this.mesh_.ghostCells_ny_bi_[i]) * this.mesh_.ghostCells_ny_bi_[i];
            this.rhorho_[ghostCell] = rho_mirror;
            this.pp_[ghostCell] = p_mirror;
            this.EE_[ghostCell] = p_mirror / ((this.inputs_.GAMMA_ - 1.0)*rho_mirror) + 0.5 * (this.uu_[ghostCell]**2 + this.vv_[ghostCell]**2);

            this.W0_[ghostCell] = this.rhorho_[ghostCell];
            this.W1_[ghostCell] = this.rhorho_[ghostCell] * this.uu_[ghostCell];
            this.W2_[ghostCell] = this.rhorho_[ghostCell] * this.vv_[ghostCell];
            this.W3_[ghostCell] = this.rhorho_[ghostCell] * this.EE_[ghostCell];

            this.ghostCells_rho_mirror_[i] = rho_mirror;
            this.ghostCells_u_mirror_[i] = u_mirror;
            this.ghostCells_v_mirror_[i] = v_mirror;
            this.ghostCells_E_mirror_[i] = E_mirror;
            this.ghostCells_p_mirror_[i] = p_mirror;
        }
        forall i in this.ghostCellIndices_dom_ {
            const ghostCell = this.ghostCellWallIndicesWithGhost_[i];
            this.ghostCells_rho_bi_[i] = 0.5 * (this.rhorho_[ghostCell] + this.ghostCells_rho_mirror_[i]);
            this.ghostCells_u_bi_[i] = 0.5 * (this.uu_[ghostCell] + this.ghostCells_u_mirror_[i]);
            this.ghostCells_v_bi_[i] = 0.5 * (this.vv_[ghostCell] + this.ghostCells_v_mirror_[i]);
            this.ghostCells_E_bi_[i] = 0.5 * (this.EE_[ghostCell] + this.ghostCells_E_mirror_[i]);
            this.ghostCells_p_bi_[i] = 0.5 * (this.pp_[ghostCell] + this.ghostCells_p_mirror_[i]);
        }
        // writeln("ghostCell; x_bi; y_bi; u_bi; v_bi; u_mirror; v_mirror; u_ghost; v_ghost");
        // for i in this.ghostCellIndices_dom_ {
        //     const ghostCell = this.ghostCellWallIndicesWithGhost_[i];
        //     const x = this.mesh_.ghostCells_x_bi_[i];
        //     const y = this.mesh_.ghostCells_y_bi_[i];
        //     writeln(ghostCell, "; ", x, "; ", y, "; ", this.ghostCells_u_bi_[i], "; ", this.ghostCells_v_bi_[i], "; ", this.ghostCells_u_mirror_[i], "; ", this.ghostCells_v_mirror_[i], "; ", this.uu_[ghostCell], "; ", this.vv_[ghostCell]);
        // }


        // Farfield
        if this.inputs_.MACH_ > 1.0 {
            forall face in this.outflowFacesI_ {
                const (i, j) = this.mesh_.IfacesIJ_[face];
                const leftCell = this.mesh_.iiCell(i-1, j);
                const leftCellWithGhost = meshIndex2FVMindex(i-1, j);
                const rightCellWithGhost = leftCellWithGhost + 1;
                this.W0_[rightCellWithGhost] = this.W0_[leftCellWithGhost];
                this.W1_[rightCellWithGhost] = this.W1_[leftCellWithGhost];
                this.W2_[rightCellWithGhost] = this.W2_[leftCellWithGhost];
                this.W3_[rightCellWithGhost] = this.W3_[leftCellWithGhost];

                this.rhorho_[rightCellWithGhost] = this.W0_[rightCellWithGhost];
                this.uu_[rightCellWithGhost] = this.W1_[rightCellWithGhost] / this.W0_[rightCellWithGhost];
                this.vv_[rightCellWithGhost] = this.W2_[rightCellWithGhost] / this.W0_[rightCellWithGhost];
                this.EE_[rightCellWithGhost] = this.W3_[rightCellWithGhost] / this.W0_[rightCellWithGhost];
                this.pp_[rightCellWithGhost] = (this.inputs_.GAMMA_ - 1.0) * this.rhorho_[rightCellWithGhost] * (this.EE_[rightCellWithGhost] - 0.5 * (this.uu_[rightCellWithGhost]**2 + this.vv_[rightCellWithGhost]**2) );
            }
            forall face in this.outflowFacesJ_ {
                const (i, j) = this.mesh_.JfacesIJ_[face];
                const bottomCell = this.mesh_.iiCell(i, j-1);
                const bottomCellWithGhost = meshIndex2FVMindex(i, j-1);
                const topCellWithGhost = bottomCellWithGhost + this.niCellWithGhosts_;
                if this.inputs_.ALPHA_ >= 0.0 {
                    this.W0_[topCellWithGhost] = this.W0_[bottomCellWithGhost];
                    this.W1_[topCellWithGhost] = this.W1_[bottomCellWithGhost];
                    this.W2_[topCellWithGhost] = this.W2_[bottomCellWithGhost];
                    this.W3_[topCellWithGhost] = this.W3_[bottomCellWithGhost];

                    this.rhorho_[topCellWithGhost] = this.W0_[topCellWithGhost];
                    this.uu_[topCellWithGhost] = this.W1_[topCellWithGhost] / this.W0_[topCellWithGhost];
                    this.vv_[topCellWithGhost] = this.W2_[topCellWithGhost] / this.W0_[topCellWithGhost];
                    this.EE_[topCellWithGhost] = this.W3_[topCellWithGhost] / this.W0_[topCellWithGhost];
                    this.pp_[topCellWithGhost] = (this.inputs_.GAMMA_ - 1.0) * this.rhorho_[topCellWithGhost] * (this.EE_[topCellWithGhost] - 0.5 * (this.uu_[topCellWithGhost]**2 + this.vv_[topCellWithGhost]**2) );
                }
                else if this.inputs_.ALPHA_ < 0.0 {
                    this.W0_[bottomCellWithGhost] = this.W0_[topCellWithGhost];
                    this.W1_[bottomCellWithGhost] = this.W1_[topCellWithGhost];
                    this.W2_[bottomCellWithGhost] = this.W2_[topCellWithGhost];
                    this.W3_[bottomCellWithGhost] = this.W3_[topCellWithGhost];

                    this.rhorho_[bottomCellWithGhost] = this.W0_[bottomCellWithGhost];
                    this.uu_[bottomCellWithGhost] = this.W1_[bottomCellWithGhost] / this.W0_[bottomCellWithGhost];
                    this.vv_[bottomCellWithGhost] = this.W2_[bottomCellWithGhost] / this.W0_[bottomCellWithGhost];
                    this.EE_[bottomCellWithGhost] = this.W3_[bottomCellWithGhost] / this.W0_[bottomCellWithGhost];
                    this.pp_[bottomCellWithGhost] = (this.inputs_.GAMMA_ - 1.0) * this.rhorho_[bottomCellWithGhost] * (this.EE_[bottomCellWithGhost] - 0.5 * (this.uu_[bottomCellWithGhost]**2 + this.vv_[bottomCellWithGhost]**2) );
                }
            }
            forall face in this.inflowFacesI_ {
                const (i, j) = this.mesh_.IfacesIJ_[face];
                const rightCellWithGhost = meshIndex2FVMindex(i, j);
                const leftCellWithGhost = rightCellWithGhost - 1;
                this.W0_[leftCellWithGhost] = 2*this.inputs_.RHO_INF_ - this.W0_[rightCellWithGhost];
                this.W1_[leftCellWithGhost] = 2*this.inputs_.RHO_INF_ * this.inputs_.U_INF_ - this.W1_[rightCellWithGhost];
                this.W2_[leftCellWithGhost] = 2*this.inputs_.RHO_INF_ * this.inputs_.V_INF_ - this.W2_[rightCellWithGhost];
                this.W3_[leftCellWithGhost] = 2*this.inputs_.RHO_INF_ * this.inputs_.E_INF_ - this.W3_[rightCellWithGhost];

                this.rhorho_[leftCellWithGhost] = this.W0_[leftCellWithGhost];
                this.uu_[leftCellWithGhost] = this.W1_[leftCellWithGhost] / this.W0_[leftCellWithGhost];
                this.vv_[leftCellWithGhost] = this.W2_[leftCellWithGhost] / this.W0_[leftCellWithGhost];
                this.EE_[leftCellWithGhost] = this.W3_[leftCellWithGhost] / this.W0_[leftCellWithGhost];
                this.pp_[leftCellWithGhost] = (this.inputs_.GAMMA_ - 1.0) * this.rhorho_[leftCellWithGhost] * (this.EE_[leftCellWithGhost] - 0.5 * (this.uu_[leftCellWithGhost]**2 + this.vv_[leftCellWithGhost]**2) );
            }
            for face in this.inflowFacesJ_ {
                const (i, j) = this.mesh_.JfacesIJ_[face];
                const topCellWithGhost = meshIndex2FVMindex(i, j);
                const bottomCellWithGhost = topCellWithGhost - this.niCellWithGhosts_;
                if this.inputs_.ALPHA_ >= 0.0 {
                    this.W0_[bottomCellWithGhost] = 2*this.inputs_.RHO_INF_ - this.W0_[topCellWithGhost];
                    this.W1_[bottomCellWithGhost] = 2*this.inputs_.RHO_INF_ * this.inputs_.U_INF_ - this.W1_[topCellWithGhost];
                    this.W2_[bottomCellWithGhost] = 2*this.inputs_.RHO_INF_ * this.inputs_.V_INF_ - this.W2_[topCellWithGhost];
                    this.W3_[bottomCellWithGhost] = 2*this.inputs_.RHO_INF_ * this.inputs_.E_INF_ - this.W3_[topCellWithGhost];

                    this.rhorho_[bottomCellWithGhost] = this.W0_[bottomCellWithGhost];
                    this.uu_[bottomCellWithGhost] = this.W1_[bottomCellWithGhost] / this.W0_[bottomCellWithGhost];
                    this.vv_[bottomCellWithGhost] = this.W2_[bottomCellWithGhost] / this.W0_[bottomCellWithGhost];
                    this.EE_[bottomCellWithGhost] = this.W3_[bottomCellWithGhost] / this.W0_[bottomCellWithGhost];
                    this.pp_[bottomCellWithGhost] = (this.inputs_.GAMMA_ - 1.0) * this.rhorho_[bottomCellWithGhost] * (this.EE_[bottomCellWithGhost] - 0.5 * (this.uu_[bottomCellWithGhost]**2 + this.vv_[bottomCellWithGhost]**2) );
                }
                else if this.inputs_.ALPHA_ < 0.0 {
                    this.W0_[topCellWithGhost] = 2*this.inputs_.RHO_INF_ - this.W0_[bottomCellWithGhost];
                    this.W1_[topCellWithGhost] = 2*this.inputs_.RHO_INF_ * this.inputs_.U_INF_ - this.W1_[bottomCellWithGhost];
                    this.W2_[topCellWithGhost] = 2*this.inputs_.RHO_INF_ * this.inputs_.V_INF_ - this.W2_[bottomCellWithGhost];
                    this.W3_[topCellWithGhost] = 2*this.inputs_.RHO_INF_ * this.inputs_.E_INF_ - this.W3_[bottomCellWithGhost];

                    this.rhorho_[topCellWithGhost] = this.W0_[topCellWithGhost];
                    this.uu_[topCellWithGhost] = this.W1_[topCellWithGhost] / this.W0_[topCellWithGhost];
                    this.vv_[topCellWithGhost] = this.W2_[topCellWithGhost] / this.W0_[topCellWithGhost];
                    this.EE_[topCellWithGhost] = this.W3_[topCellWithGhost] / this.W0_[topCellWithGhost];
                    this.pp_[topCellWithGhost] = (this.inputs_.GAMMA_ - 1.0) * this.rhorho_[topCellWithGhost] * (this.EE_[topCellWithGhost] - 0.5 * (this.uu_[topCellWithGhost]**2 + this.vv_[topCellWithGhost]**2) );
                }
            }
        }
    }

    proc updatePrimitiveVariablesFromConserved() {
        forall cell in this.cell_dom_with_ghosts_ {
            this.rhorho_[cell] = this.W0_[cell];
            this.uu_[cell] = this.W1_[cell] / this.W0_[cell];
            this.vv_[cell] = this.W2_[cell] / this.W0_[cell];
            this.EE_[cell] = this.W3_[cell] / this.W0_[cell];
            this.pp_[cell] = (this.inputs_.GAMMA_ - 1.0) * this.rhorho_[cell] * (this.EE_[cell] - 0.5 * (this.uu_[cell]**2 + this.vv_[cell]**2) );
        }
    }

    proc compute_convective_fluxes() {
        forall face in this.mesh_.face_dom_ {
            const (i, j) = this.mesh_.IfacesIJ_[face];
            const leftCell = meshIndex2FVMindex(i-1, j);
            const rightCell = meshIndex2FVMindex(i, j);

            if (this.cellTypesWithGhosts_[leftCell] == 0 && this.cellTypesWithGhosts_[rightCell] == 0) {
                continue;
            }

            const avg_W0 = 0.5 * (this.W0_[leftCell] + this.W0_[rightCell]);
            const avg_W1 = 0.5 * (this.W1_[leftCell] + this.W1_[rightCell]);
            const avg_W2 = 0.5 * (this.W2_[leftCell] + this.W2_[rightCell]);
            const avg_W3 = 0.5 * (this.W3_[leftCell] + this.W3_[rightCell]);

            const avg_u = 0.5 * (this.uu_[leftCell] + this.uu_[rightCell]);
            const avg_v = 0.5 * (this.vv_[leftCell] + this.vv_[rightCell]);
            const avg_E = 0.5 * (this.EE_[leftCell] + this.EE_[rightCell]);
            const avg_p = 0.5 * (this.pp_[leftCell] + this.pp_[rightCell]);

            this.F0I_[face] = avg_W0 * avg_u * this.mesh_.IfaceAreas_[face];
            this.F1I_[face] = (avg_W1 * avg_u + avg_p) * this.mesh_.IfaceAreas_[face];
            this.F2I_[face] = (avg_W2 * avg_u) * this.mesh_.IfaceAreas_[face];
            this.F3I_[face] = (avg_W3 * avg_u + avg_p * avg_u) * this.mesh_.IfaceAreas_[face];
        }

        forall face in this.mesh_.face_dom_ {
            const (i, j) = this.mesh_.JfacesIJ_[face];
            const bottomCell = meshIndex2FVMindex(i, j-1);
            const topCell = meshIndex2FVMindex(i, j);

            if (this.cellTypesWithGhosts_[bottomCell] == 0 && this.cellTypesWithGhosts_[topCell] == 0) {
                continue;
            }

            const avg_W0 = 0.5 * (this.W0_[bottomCell] + this.W0_[topCell]);
            const avg_W1 = 0.5 * (this.W1_[bottomCell] + this.W1_[topCell]);
            const avg_W2 = 0.5 * (this.W2_[bottomCell] + this.W2_[topCell]);
            const avg_W3 = 0.5 * (this.W3_[bottomCell] + this.W3_[topCell]);

            const avg_u = 0.5 * (this.uu_[bottomCell] + this.uu_[topCell]);
            const avg_v = 0.5 * (this.vv_[bottomCell] + this.vv_[topCell]);
            const avg_E = 0.5 * (this.EE_[bottomCell] + this.EE_[topCell]);
            const avg_p = 0.5 * (this.pp_[bottomCell] + this.pp_[topCell]);

            this.F0J_[face] = avg_W0 * avg_v * this.mesh_.JfaceAreas_[face];
            this.F1J_[face] = (avg_W1 * avg_v) * this.mesh_.JfaceAreas_[face];
            this.F2J_[face] = (avg_W2 * avg_v + avg_p) * this.mesh_.JfaceAreas_[face];
            this.F3J_[face] = (avg_W3 * avg_v + avg_p * avg_v) * this.mesh_.JfaceAreas_[face];
        }
    }

    proc compute_lambdas() {
        forall j in 0..<this.mesh_.njCell_ {
            for i in 0..<this.mesh_.niCell_ {
                const idx = this.mesh_.iiCell(i, j);
                if this.mesh_.cellTypes_[idx] == 0 {
                    continue;
                }
                const idxWithGhost = meshIndex2FVMindex(i, j);
                const u = this.uu_[idxWithGhost];
                const v = this.vv_[idxWithGhost];
                const p = this.pp_[idxWithGhost];
                const rho = this.rhorho_[idxWithGhost];
                const a = sqrt(this.inputs_.GAMMA_ * p / rho);

                this.LambdaI_[idxWithGhost] = (abs(u) + a) * this.mesh_.avgFaceAreaI_[idx];
                this.LambdaJ_[idxWithGhost] = (abs(v) + a) * this.mesh_.avgFaceAreaJ_[idx];
            }
        }
    }

    proc meshIndex2FVMindex(i: int, j: int): int {
        return i+2 + (j+2) * this.niCellWithGhosts_;
    }
}





}