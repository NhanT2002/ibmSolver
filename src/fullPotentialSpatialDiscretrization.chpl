module fullPotentialSpatialDiscretization 
{
use mesh;
use writeCGNS;
use Time;
use Random;
use List;
import input.inputsConfig;
use Math;
use Sort;
use Map;

class fullPotentialSpatialDiscretization {
    var time_: stopwatch;
    var mesh_ : shared meshData;
    var inputs_ : inputsConfig;

    var one_over_gamma_ : real(64);
    var one_over_gamma_minus_one_ : real(64);
    var two_gamma_over_gamma_minus_one_ : real(64);

    var cell_dom_with_ghosts_ : domain(1) = {1..0};
    var nCellWithGhosts_ : int;
    var niCellWithGhosts_ : int;
    var njCellWithGhosts_ : int;
    var xCellsWithGhosts_ : [cell_dom_with_ghosts_] real(64);
    var yCellsWithGhosts_ : [cell_dom_with_ghosts_] real(64);
    var cellTypesWithGhosts_ : [cell_dom_with_ghosts_] int;
    var cellVolumesWithGhosts_ : [cell_dom_with_ghosts_] real(64);

    var rhorho_ : [cell_dom_with_ghosts_] real(64);
    var phiphi_ : [cell_dom_with_ghosts_] real(64);
    var uu_ : [cell_dom_with_ghosts_] real(64);
    var vv_ : [cell_dom_with_ghosts_] real(64);
    var pp_ : [cell_dom_with_ghosts_] real(64);

    var mi_ : [cell_dom_with_ghosts_] real(64);
    var si_ : [cell_dom_with_ghosts_] real(64);
    var ni_ : [cell_dom_with_ghosts_] real(64);

    var mj_ : [cell_dom_with_ghosts_] real(64);
    var sj_ : [cell_dom_with_ghosts_] real(64);
    var nj_ : [cell_dom_with_ghosts_] real(64);

    var Rc0_ : [cell_dom_with_ghosts_] real(64);
    var Rc1_ : [cell_dom_with_ghosts_] real(64);
    var phi_t_ : [cell_dom_with_ghosts_] real(64);

    var LambdaI_ : [cell_dom_with_ghosts_] real(64);
    var LambdaJ_ : [cell_dom_with_ghosts_] real(64);

    var Rd0_ : [cell_dom_with_ghosts_] real(64);
    var Rd1_ : [cell_dom_with_ghosts_] real(64);
    var Rd2_ : [cell_dom_with_ghosts_] real(64);
    var Rd3_ : [cell_dom_with_ghosts_] real(64);

    var R0_ : [cell_dom_with_ghosts_] real(64);
    var R1_ : [cell_dom_with_ghosts_] real(64);
    var R2_ : [cell_dom_with_ghosts_] real(64);
    var R3_ : [cell_dom_with_ghosts_] real(64);

    // 1st layer of ghost cells
    var ghostCellIndices_dom_ : domain(1) = {1..0};
    var ghostCellWallIndicesWithGhost_ : [ghostCellIndices_dom_] int;
    var ghostCellsNearestFluidCellsWithGhost_dom_ : domain(2) = {1..0, 1..0};
    var ghostCellsNearestFluidCellsWithGhost_ : [ghostCellsNearestFluidCellsWithGhost_dom_] int;

    var ghostCells_rho_bi_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_phi_bi_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_u_bi_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_v_bi_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_w_bi_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_cp_bi_ : [ghostCellIndices_dom_] real(64);

    var ghostCells_rho_mirror_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_phi_mirror_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_u_mirror_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_v_mirror_ : [ghostCellIndices_dom_] real(64);
    var ghostCells_p_mirror_ : [ghostCellIndices_dom_] real(64);

    // 2nd layer of ghost cells
    var ghostCellm1Indices_dom_ : domain(1) = {1..0};
    var ghostCellm1WallIndicesWithGhost_ : [ghostCellm1Indices_dom_] int;
    var ghostCellsm1NearestFluidCellsWithGhost_dom_ : domain(2) = {1..0, 1..0};
    var ghostCellsm1NearestFluidCellsWithGhost_ : [ghostCellsm1NearestFluidCellsWithGhost_dom_] int;

    var ghostCellsm1_rho_bi_ : [ghostCellm1Indices_dom_] real(64);
    var ghostCellsm1_phi_bi_ : [ghostCellm1Indices_dom_] real(64);
    var ghostCellsm1_u_bi_ : [ghostCellm1Indices_dom_] real(64);
    var ghostCellsm1_v_bi_ : [ghostCellm1Indices_dom_] real(64);
    var ghostCellsm1_w_bi_ : [ghostCellm1Indices_dom_] real(64);
    var ghostCellsm1_cp_bi_ : [ghostCellm1Indices_dom_] real(64);

    var ghostCellsm1_rho_mirror_ : [ghostCellm1Indices_dom_] real(64);
    var ghostCellsm1_phi_mirror_ : [ghostCellm1Indices_dom_] real(64);
    var ghostCellsm1_u_mirror_ : [ghostCellm1Indices_dom_] real(64);
    var ghostCellsm1_v_mirror_ : [ghostCellm1Indices_dom_] real(64);
    var ghostCellsm1_p_mirror_ : [ghostCellm1Indices_dom_] real(64);

    var face_dom_ : domain(1) = {1..0};
    var F0I_ : [face_dom_] real(64);
    var F1I_ : [face_dom_] real(64);
    var dF0Idrho_ : [face_dom_] real(64);
    var avg_rho_I_ : [face_dom_] real(64);

    var F0J_ : [face_dom_] real(64);
    var F1J_ : [face_dom_] real(64);
    var dF0Jdrho_ : [face_dom_] real(64);
    var avg_rho_J_ : [face_dom_] real(64);

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

        this.one_over_gamma_ = 1.0 / this.inputs_.GAMMA_;
        this.one_over_gamma_minus_one_ = 1.0 / (this.inputs_.GAMMA_ - 1.0);
        this.two_gamma_over_gamma_minus_one_ = 2.0 * this.inputs_.GAMMA_ / (this.inputs_.GAMMA_ - 1.0);

        // Define domains
        const niCell_ = mesh_.niCell_;
        const njCell_ = mesh_.njCell_;
        this.cell_dom_with_ghosts_ = {0..<( (niCell_ + 4) * (njCell_ + 4) )};
        this.nCellWithGhosts_ = (niCell_ + 4) * (njCell_ + 4);
        this.niCellWithGhosts_ = niCell_ + 4;
        this.njCellWithGhosts_ = njCell_ + 4;
        this.ghostCellIndices_dom_ = this.mesh_.ghostCellIndices_.domain;
        this.ghostCellsNearestFluidCellsWithGhost_dom_ = this.mesh_.ghostCellsNearestFluidCells_.domain;
        this.ghostCellm1Indices_dom_ = this.mesh_.ghostCellm1Indices_.domain;
        this.ghostCellsm1NearestFluidCellsWithGhost_dom_ = this.mesh_.ghostCellsm1NearestFluidCells_.domain;
        this.face_dom_ = {0..<( (niCell_+1)*njCell_ )};

        writeln("spatialDiscretization initialized.");
        writeln("  niCell_ = ", niCell_);
        writeln("  njCell_ = ", njCell_);
        writeln("  niCellWithGhosts_ = ", this.niCellWithGhosts_);
        writeln("  njCellWithGhosts_ = ", this.njCellWithGhosts_);
        writeln("  MACH_ = ", this.inputs_.MACH_);
        writeln("  ALPHA_ = ", this.inputs_.ALPHA_);
        writeln("  U_INF_ = ", this.inputs_.U_INF_);
        writeln("  V_INF_ = ", this.inputs_.V_INF_);
    }

    proc initializeFlowField() {

        // Map to mesh cells with ghosts
        this.xCellsWithGhosts_ = 0;
        this.yCellsWithGhosts_ = 0;
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
            const face = i-2;
            this.xCellsWithGhosts_[i+this.niCellWithGhosts_] = this.xCellsWithGhosts_[i + 2 * this.niCellWithGhosts_];
            this.yCellsWithGhosts_[i+this.niCellWithGhosts_] = this.yCellsWithGhosts_[i + 2 * this.niCellWithGhosts_] - 2*abs(this.yCellsWithGhosts_[i + 2 * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[face]);
        }
        forall i in 2..<this.niCellWithGhosts_-2 {
            const face = i-2;
            this.xCellsWithGhosts_[i] = this.xCellsWithGhosts_[i+this.niCellWithGhosts_];
            this.yCellsWithGhosts_[i] = this.yCellsWithGhosts_[i+this.niCellWithGhosts_] - 2*abs(this.yCellsWithGhosts_[i + 2 * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[face]);
        }
        forall i in 2..<this.niCellWithGhosts_-2 {
            const face = i-2 + (this.mesh_.njCell_)* (this.mesh_.niCell_);
            this.xCellsWithGhosts_[i + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] = this.xCellsWithGhosts_[i + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_];
            this.yCellsWithGhosts_[i + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] = this.yCellsWithGhosts_[i + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_] + 2*abs(this.yCellsWithGhosts_[i + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[face]);
        }
        forall i in 2..<this.niCellWithGhosts_-2 {
            const face = i-2 + (this.mesh_.njCell_)* (this.mesh_.niCell_);
            this.xCellsWithGhosts_[i + (this.njCellWithGhosts_-1) * this.niCellWithGhosts_] = this.xCellsWithGhosts_[i + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_];
            this.yCellsWithGhosts_[i + (this.njCellWithGhosts_-1) * this.niCellWithGhosts_] = this.yCellsWithGhosts_[i + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] + 2*abs(this.yCellsWithGhosts_[i + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[face]);
        }
        forall j in 2..<this.njCellWithGhosts_-2 {
            const face = (j-2) * (this.mesh_.niCell_+1);
            this.xCellsWithGhosts_[1 + j * this.niCellWithGhosts_] = this.xCellsWithGhosts_[2 + j * this.niCellWithGhosts_] - 2*abs(this.xCellsWithGhosts_[2 + j * this.niCellWithGhosts_] - this.mesh_.IfacesCx_[face]);
            this.yCellsWithGhosts_[1 + j * this.niCellWithGhosts_] = this.yCellsWithGhosts_[2 + j * this.niCellWithGhosts_];
        }
        forall j in 2..<this.njCellWithGhosts_-2 {
            const face = (j-2) * (this.mesh_.niCell_+1);
            this.xCellsWithGhosts_[0 + j * this.niCellWithGhosts_] = this.xCellsWithGhosts_[1 + j * this.niCellWithGhosts_] - 2*abs(this.xCellsWithGhosts_[2 + j * this.niCellWithGhosts_] - this.mesh_.IfacesCx_[face]);
            this.yCellsWithGhosts_[0 + j * this.niCellWithGhosts_] = this.yCellsWithGhosts_[1 + j * this.niCellWithGhosts_];
        }
        forall j in 2..<this.njCellWithGhosts_-2 {
            const face = this.mesh_.niCell_ + (j-2) * (this.mesh_.niCell_+1);
            this.xCellsWithGhosts_[(this.niCellWithGhosts_-2) + j * this.niCellWithGhosts_] = this.xCellsWithGhosts_[(this.niCellWithGhosts_-3) + j * this.niCellWithGhosts_] + 2*abs(this.xCellsWithGhosts_[(this.niCellWithGhosts_-3) + j * this.niCellWithGhosts_] - this.mesh_.IfacesCx_[face]);
            this.yCellsWithGhosts_[(this.niCellWithGhosts_-2) + j * this.niCellWithGhosts_] = this.yCellsWithGhosts_[(this.niCellWithGhosts_-3) + j * this.niCellWithGhosts_];
        }
        forall j in 2..<this.njCellWithGhosts_-2 {
            const face = this.mesh_.niCell_ + (j-2) * (this.mesh_.niCell_+1);
            this.xCellsWithGhosts_[(this.niCellWithGhosts_-1) + j * this.niCellWithGhosts_] = this.xCellsWithGhosts_[(this.niCellWithGhosts_-2) + j * this.niCellWithGhosts_] + 2*abs(this.xCellsWithGhosts_[(this.niCellWithGhosts_-3) + j * this.niCellWithGhosts_] - this.mesh_.IfacesCx_[face]);
            this.yCellsWithGhosts_[(this.niCellWithGhosts_-1) + j * this.niCellWithGhosts_] = this.yCellsWithGhosts_[(this.niCellWithGhosts_-2) + j * this.niCellWithGhosts_];
        }
        // bottom left corner cells
        this.xCellsWithGhosts_[0 + 1 * this.niCellWithGhosts_] = this.xCellsWithGhosts_[0 + 2 * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[0 + 1 * this.niCellWithGhosts_] = this.yCellsWithGhosts_[0 + 2 * this.niCellWithGhosts_] - 2*abs(this.yCellsWithGhosts_[0 + 2 * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[0]);
        this.xCellsWithGhosts_[1 + 1 * this.niCellWithGhosts_] = this.xCellsWithGhosts_[1 + 2 * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[1 + 1 * this.niCellWithGhosts_] = this.yCellsWithGhosts_[1 + 2 * this.niCellWithGhosts_] - 2*abs(this.yCellsWithGhosts_[1 + 2 * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[0]);

        this.xCellsWithGhosts_[0 + 0 * this.niCellWithGhosts_] = this.xCellsWithGhosts_[0 + 1 * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[0 + 0 * this.niCellWithGhosts_] = this.yCellsWithGhosts_[0 + 1 * this.niCellWithGhosts_] - 2*abs(this.yCellsWithGhosts_[0 + 2 * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[0]);
        this.xCellsWithGhosts_[1 + 0 * this.niCellWithGhosts_] = this.xCellsWithGhosts_[1 + 1 * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[1 + 0 * this.niCellWithGhosts_] = this.yCellsWithGhosts_[1 + 1 * this.niCellWithGhosts_] - 2*abs(this.yCellsWithGhosts_[1 + 2 * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[0]);
        
        // bottom right corner cells
        this.xCellsWithGhosts_[(this.niCellWithGhosts_-2) + 1 * this.niCellWithGhosts_] = this.xCellsWithGhosts_[(this.niCellWithGhosts_-2) + 2 * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[(this.niCellWithGhosts_-2) + 1 * this.niCellWithGhosts_] = this.yCellsWithGhosts_[(this.niCellWithGhosts_-2) + 2 * this.niCellWithGhosts_] - 2*abs(this.yCellsWithGhosts_[(this.niCellWithGhosts_-2) + 2 * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[(this.niCellWithGhosts_-5)]); 
        this.xCellsWithGhosts_[(this.niCellWithGhosts_-1) + 1 * this.niCellWithGhosts_] = this.xCellsWithGhosts_[(this.niCellWithGhosts_-1) + 2 * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[(this.niCellWithGhosts_-1) + 1 * this.niCellWithGhosts_] = this.yCellsWithGhosts_[(this.niCellWithGhosts_-1) + 2 * this.niCellWithGhosts_] - 2*abs(this.yCellsWithGhosts_[(this.niCellWithGhosts_-1) + 2 * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[(this.niCellWithGhosts_-5)]);

        this.xCellsWithGhosts_[(this.niCellWithGhosts_-2) + 0 * this.niCellWithGhosts_] = this.xCellsWithGhosts_[(this.niCellWithGhosts_-2) + 1 * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[(this.niCellWithGhosts_-2) + 0 * this.niCellWithGhosts_] = this.yCellsWithGhosts_[(this.niCellWithGhosts_-2) + 1 * this.niCellWithGhosts_] - 2*abs(this.yCellsWithGhosts_[(this.niCellWithGhosts_-2) + 2 * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[(this.niCellWithGhosts_-5)]);
        this.xCellsWithGhosts_[(this.niCellWithGhosts_-1) + 0 * this.niCellWithGhosts_] = this.xCellsWithGhosts_[(this.niCellWithGhosts_-1) + 1 * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[(this.niCellWithGhosts_-1) + 0 * this.niCellWithGhosts_] = this.yCellsWithGhosts_[(this.niCellWithGhosts_-1) + 1 * this.niCellWithGhosts_] - 2*abs(this.yCellsWithGhosts_[(this.niCellWithGhosts_-1) + 2 * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[(this.niCellWithGhosts_-5)]);

        // top left corner cells
        this.xCellsWithGhosts_[0 + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] = this.xCellsWithGhosts_[0 + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[0 + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] = this.yCellsWithGhosts_[0 + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_] + 2*abs(this.yCellsWithGhosts_[0 + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[0 + (this.njCellWithGhosts_-5) * (this.mesh_.niCell_)]);
        this.xCellsWithGhosts_[1 + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] = this.xCellsWithGhosts_[1 + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[1 + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] = this.yCellsWithGhosts_[1 + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_] + 2*abs(this.yCellsWithGhosts_[1 + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[0 + (this.njCellWithGhosts_-5) * (this.mesh_.niCell_)]);

        this.xCellsWithGhosts_[0 + (this.njCellWithGhosts_-1) * this.niCellWithGhosts_] = this.xCellsWithGhosts_[0 + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[0 + (this.njCellWithGhosts_-1) * this.niCellWithGhosts_] = this.yCellsWithGhosts_[0 + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] + 2*abs(this.yCellsWithGhosts_[0 + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[0 + (this.njCellWithGhosts_-5) * (this.mesh_.niCell_)]);
        this.xCellsWithGhosts_[1 + (this.njCellWithGhosts_-1) * this.niCellWithGhosts_] = this.xCellsWithGhosts_[1 + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[1 + (this.njCellWithGhosts_-1) * this.niCellWithGhosts_] = this.yCellsWithGhosts_[1 + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] + 2*abs(this.yCellsWithGhosts_[1 + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[0 + (this.njCellWithGhosts_-5) * (this.mesh_.niCell_)]);

        // top right corner cells
        this.xCellsWithGhosts_[(this.niCellWithGhosts_-2) + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] = this.xCellsWithGhosts_[(this.niCellWithGhosts_-2) + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[(this.niCellWithGhosts_-2) + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] = this.yCellsWithGhosts_[(this.niCellWithGhosts_-2) + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_] + 2*abs(this.yCellsWithGhosts_[(this.niCellWithGhosts_-2) + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[(this.niCellWithGhosts_-5) + (this.njCellWithGhosts_-5) * (this.mesh_.niCell_)]);
        this.xCellsWithGhosts_[(this.niCellWithGhosts_-1) + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] = this.xCellsWithGhosts_[(this.niCellWithGhosts_-1) + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[(this.niCellWithGhosts_-1) + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] = this.yCellsWithGhosts_[(this.niCellWithGhosts_-1) + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_] + 2*abs(this.yCellsWithGhosts_[(this.niCellWithGhosts_-1) + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[(this.niCellWithGhosts_-5) + (this.njCellWithGhosts_-5) * (this.mesh_.niCell_)]);

        this.xCellsWithGhosts_[(this.niCellWithGhosts_-2) + (this.njCellWithGhosts_-1) * this.niCellWithGhosts_] = this.xCellsWithGhosts_[(this.niCellWithGhosts_-2) + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[(this.niCellWithGhosts_-2) + (this.njCellWithGhosts_-1) * this.niCellWithGhosts_] = this.yCellsWithGhosts_[(this.niCellWithGhosts_-2) + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] + 2*abs(this.yCellsWithGhosts_[(this.niCellWithGhosts_-2) + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[(this.niCellWithGhosts_-5) + (this.njCellWithGhosts_-5) * (this.mesh_.niCell_)]);
        this.xCellsWithGhosts_[(this.niCellWithGhosts_-1) + (this.njCellWithGhosts_-1) * this.niCellWithGhosts_] = this.xCellsWithGhosts_[(this.niCellWithGhosts_-1) + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_];
        this.yCellsWithGhosts_[(this.niCellWithGhosts_-1) + (this.njCellWithGhosts_-1) * this.niCellWithGhosts_] = this.yCellsWithGhosts_[(this.niCellWithGhosts_-1) + (this.njCellWithGhosts_-2) * this.niCellWithGhosts_] + 2*abs(this.yCellsWithGhosts_[(this.niCellWithGhosts_-1) + (this.njCellWithGhosts_-3) * this.niCellWithGhosts_] - this.mesh_.JfacesCy_[(this.niCellWithGhosts_-5) + (this.njCellWithGhosts_-5) * (this.mesh_.niCell_)]);

        forall (id,(i, j)) in zip(this.mesh_.ghostCellIJ_.domain, this.mesh_.ghostCellIJ_) {
            const idxWithGhost = i+2 + (j+2) * this.niCellWithGhosts_;
            this.ghostCellWallIndicesWithGhost_[id] = idxWithGhost;
        }
        forall i in this.ghostCellIndices_dom_ {
            for j in this.mesh_.ghostCellsNearestFluidCellsCx_dom_.dim(1) {
                const (fi, fj) = this.mesh_.ghostCellsNearestFluidCellsIJ_[i, j];
                const fi_withGhost = fi + 2;
                const fj_withGhost = fj + 2;
                const idxWithGhost = fi_withGhost + fj_withGhost * this.niCellWithGhosts_;
                this.ghostCellsNearestFluidCellsWithGhost_[i, j] = idxWithGhost;
            }
        }
        forall (id,(i, j)) in zip(this.mesh_.ghostCellm1IJ_.domain, this.mesh_.ghostCellm1IJ_) {
            const idxWithGhost = i+2 + (j+2) * this.niCellWithGhosts_;
            this.ghostCellm1WallIndicesWithGhost_[id] = idxWithGhost;
        }
        forall i in this.ghostCellm1Indices_dom_ {
            for j in this.mesh_.ghostCellsm1NearestFluidCellsCx_dom_.dim(1) {
                const (fi, fj) = this.mesh_.ghostCellsm1NearestFluidCellsIJ_[i, j];
                const fi_withGhost = fi + 2;
                const fj_withGhost = fj + 2;
                const idxWithGhost = fi_withGhost + fj_withGhost * this.niCellWithGhosts_;
                this.ghostCellsm1NearestFluidCellsWithGhost_[i, j] = idxWithGhost;
            }
        }

        for j in this.mesh_.nearestCellLEIJ_.domain {
            const (i_cell, j_cell) = this.mesh_.nearestCellLEIJ_[j];
            const i_cell_withGhost = i_cell + 2;
            const j_cell_withGhost = j_cell + 2;
            const idxWithGhost = i_cell_withGhost + j_cell_withGhost * this.niCellWithGhosts_;
            this.mesh_.nearestCellLE_[j] = idxWithGhost;
        }
        for j in this.mesh_.nearestCellTEIJ_.domain {
            const (i_cell, j_cell) = this.mesh_.nearestCellTEIJ_[j];
            const i_cell_withGhost = i_cell + 2;
            const j_cell_withGhost = j_cell + 2;
            const idxWithGhost = i_cell_withGhost + j_cell_withGhost * this.niCellWithGhosts_;
            this.mesh_.nearestCellTE_[j] = idxWithGhost;
        }


        this.cellTypesWithGhosts_ = 9;
        forall j in 0..<this.mesh_.njCell_ {
            for i in 0..<this.mesh_.niCell_ {
                const idxWithGhost = meshIndex2FVMindex(i, j);
                const idx = this.mesh_.iiCell(i, j);
                this.cellTypesWithGhosts_[idxWithGhost] = this.mesh_.cellTypes_[idx];
                this.cellVolumesWithGhosts_[idxWithGhost] = this.mesh_.cellVolumes_[idx];
            }
        }

        if this.inputs_.INITIAL_SOLUTION_ != "" {
            writeln("Reading initial solution from ", this.inputs_.INITIAL_SOLUTION_);
            const (Density, Phi, VelocityX, VelocityY) = readFullPotentialCGNSFlowField(this.inputs_.INITIAL_SOLUTION_);

            if Density.size != this.mesh_.nCell_ {
                halt("Error: Initial solution size does not match number of cells with ghosts.");
            }
            else {
                forall j in 0..<this.mesh_.njCell_ {
                    for i in 0..<this.mesh_.niCell_ {
                        const idxWithGhost = meshIndex2FVMindex(i, j);
                
                        this.rhorho_[idxWithGhost] = Density[j, i];
                        this.phiphi_[idxWithGhost] = Phi[j, i];
                        this.uu_[idxWithGhost] = VelocityX[j, i];
                        this.vv_[idxWithGhost] = VelocityY[j, i];
                    }
                }
            }
        }

        else {
            forall j in 0..<this.mesh_.njCell_ {
                for i in 0..<this.mesh_.niCell_ {
                    const idxWithGhost = meshIndex2FVMindex(i, j);
                    
                    this.rhorho_[idxWithGhost] = this.inputs_.RHO_INF_;
                    this.phiphi_[idxWithGhost] = this.inputs_.U_INF_ * this.xCellsWithGhosts_[idxWithGhost] + this.inputs_.V_INF_ * this.yCellsWithGhosts_[idxWithGhost];
                    this.uu_[idxWithGhost] = this.inputs_.U_INF_;
                    this.vv_[idxWithGhost] = this.inputs_.V_INF_;
                }
            }
        }

        

        this.ghostCells_rho_bi_ = this.inputs_.RHO_INF_;
        this.ghostCells_u_bi_ = this.inputs_.U_INF_;
        this.ghostCells_v_bi_ = this.inputs_.V_INF_;
        this.ghostCells_cp_bi_ = 0.0;

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

        forall j in 1..<this.njCellWithGhosts_-1 {
            for i in 1..<this.niCellWithGhosts_-1 {
                const cell = i + j * this.niCellWithGhosts_;
                const leftCell = cell - 1;
                const rightCell = cell + 1;
                const bottomCell = cell - this.niCellWithGhosts_;
                const topCell = cell + this.niCellWithGhosts_;

                const h_im1 = this.xCellsWithGhosts_[cell] - this.xCellsWithGhosts_[leftCell];
                const h_ip1 = this.xCellsWithGhosts_[rightCell] - this.xCellsWithGhosts_[cell];
                const h_jm1 = this.yCellsWithGhosts_[cell] - this.yCellsWithGhosts_[bottomCell];
                const h_jp1 = this.yCellsWithGhosts_[topCell] - this.yCellsWithGhosts_[cell];

                const m_i = -h_ip1 / (h_im1*(h_ip1 + h_im1));
                const n_i = h_im1 / (h_ip1*(h_ip1 + h_im1));
                const s_i = (h_ip1**2 - h_im1**2) / (h_ip1 * h_im1 * (h_ip1 + h_im1));

                const m_j = -h_jp1 / (h_jm1*(h_jp1 + h_jm1));
                const n_j = h_jm1 / (h_jp1*(h_jp1 + h_jm1));
                const s_j = (h_jp1**2 - h_jm1**2) / (h_jp1 * h_jm1 * (h_jp1 + h_jm1));

                this.mi_[cell] = m_i;
                this.ni_[cell] = n_i;
                this.si_[cell] = s_i;

                this.mj_[cell] = m_j;
                this.nj_[cell] = n_j;
                this.sj_[cell] = s_j;
            }
        }

        // forall cell in 0..<this.nCellWithGhosts_{
        //     if this.cellTypesWithGhosts_[cell] != 1 {
        //         continue;
        //     }
        //     const leftCell = cell - 1;
        //     const rightCell = cell + 1;
        //     const bottomCell = cell - this.niCellWithGhosts_;
        //     const topCell = cell + this.niCellWithGhosts_;

        //     const h_im1 = this.xCellsWithGhosts_[leftCell] - this.xCellsWithGhosts_[cell];
        //     const h_ip1 = this.xCellsWithGhosts_[rightCell] - this.xCellsWithGhosts_[cell];
        //     const h_jm1 = this.yCellsWithGhosts_[bottomCell] - this.yCellsWithGhosts_[cell];
        //     const h_jp1 = this.yCellsWithGhosts_[topCell] - this.yCellsWithGhosts_[cell];

        //     const m_i = -h_ip1 / (h_im1**2 - h_ip1*h_im1);
        //     const n_i = -h_im1 / (h_ip1**2 - h_ip1*h_im1);
        //     const s_i = -(m_i + n_i);

        //     const m_j = -h_jp1 / (h_jm1**2 - h_jp1*h_jm1);
        //     const n_j = -h_jm1 / (h_jp1**2 - h_jp1*h_jm1);
        //     const s_j = -(m_j + n_j);

        //     this.mi_[cell] = m_i;
        //     this.ni_[cell] = n_i;
        //     this.si_[cell] = s_i;

        //     this.mj_[cell] = m_j;
        //     this.nj_[cell] = n_j;
        //     this.sj_[cell] = s_j;
        // }
    }

    proc updateGhostCells() {
        // Wall
        forall i in 0..<this.ghostCellWallIndicesWithGhost_.size {
            const ghostCell = this.ghostCellWallIndicesWithGhost_[i];
            var kls = this.mesh_.ghostCellkls_[i]!.borrow();
            for (j, nearestCell) in zip(0..<4, this.ghostCellsNearestFluidCellsWithGhost_[i, 0..]) { // MAYBE GENERALIZE LATER
                kls.rho_fieldValues_[j] = this.rhorho_[nearestCell];
                kls.phi_fieldValues_[j] = this.phiphi_[nearestCell];
                kls.u_fieldValues_[j] = this.uu_[nearestCell];
                kls.v_fieldValues_[j] = this.vv_[nearestCell];
            }

            const rho_mirror = kls.interpolate(kls.rho_fieldValues_);
            const phi_mirror = kls.interpolate(kls.phi_fieldValues_);
            const u_mirror = kls.interpolate(kls.u_fieldValues_);
            const v_mirror = kls.interpolate(kls.v_fieldValues_);

            this.rhorho_[ghostCell] = rho_mirror;
            this.phiphi_[ghostCell] = phi_mirror;
            this.uu_[ghostCell] = u_mirror - 2.0 * (u_mirror * this.mesh_.ghostCells_nx_bi_[i] + v_mirror * this.mesh_.ghostCells_ny_bi_[i]) * this.mesh_.ghostCells_nx_bi_[i];
            this.vv_[ghostCell] = v_mirror - 2.0 * (u_mirror * this.mesh_.ghostCells_nx_bi_[i] + v_mirror * this.mesh_.ghostCells_ny_bi_[i]) * this.mesh_.ghostCells_ny_bi_[i];

            this.ghostCells_rho_mirror_[i] = rho_mirror;
            this.ghostCells_phi_mirror_[i] = phi_mirror;
            this.ghostCells_u_mirror_[i] = u_mirror;
            this.ghostCells_v_mirror_[i] = v_mirror;
        }
        forall i in 0..<this.ghostCellWallIndicesWithGhost_.size {
            const ghostCell = this.ghostCellWallIndicesWithGhost_[i];
            this.ghostCells_rho_bi_[i] = 0.5 * (this.rhorho_[ghostCell] + this.ghostCells_rho_mirror_[i]);
            this.ghostCells_phi_bi_[i] = 0.5 * (this.phiphi_[ghostCell] + this.ghostCells_phi_mirror_[i]);
            this.ghostCells_u_bi_[i] = 0.5 * (this.uu_[ghostCell] + this.ghostCells_u_mirror_[i]);
            this.ghostCells_v_bi_[i] = 0.5 * (this.vv_[ghostCell] + this.ghostCells_v_mirror_[i]);
            this.ghostCells_cp_bi_[i] = 2*(this.ghostCells_rho_bi_[i]**this.inputs_.GAMMA_ - 1.0) / (this.inputs_.GAMMA_ * this.inputs_.MACH_**2);
        }

        // wall m-1
        forall i in 0..<this.ghostCellm1WallIndicesWithGhost_.size {
            const ghostCell = this.ghostCellm1WallIndicesWithGhost_[i];
            var kls = this.mesh_.ghostCellsm1kls_[i]!.borrow();
            for (j, nearestCell) in zip(0..<4, this.ghostCellsm1NearestFluidCellsWithGhost_[i, 0..]) { // MAYBE GENERALIZE LATER
                kls.rho_fieldValues_[j] = this.rhorho_[nearestCell];
                kls.phi_fieldValues_[j] = this.phiphi_[nearestCell];
                kls.u_fieldValues_[j] = this.uu_[nearestCell];
                kls.v_fieldValues_[j] = this.vv_[nearestCell];
            }

            const rho_mirror = kls.interpolate(kls.rho_fieldValues_);
            const phi_mirror = kls.interpolate(kls.phi_fieldValues_);
            const u_mirror = kls.interpolate(kls.u_fieldValues_);
            const v_mirror = kls.interpolate(kls.v_fieldValues_);

            this.rhorho_[ghostCell] = rho_mirror;
            this.phiphi_[ghostCell] = phi_mirror;
            this.uu_[ghostCell] = u_mirror - 2.0 * (u_mirror * this.mesh_.ghostCellsm1_nx_bi_[i] + v_mirror * this.mesh_.ghostCellsm1_ny_bi_[i]) * this.mesh_.ghostCellsm1_nx_bi_[i];
            this.vv_[ghostCell] = v_mirror - 2.0 * (u_mirror * this.mesh_.ghostCellsm1_nx_bi_[i] + v_mirror * this.mesh_.ghostCellsm1_ny_bi_[i]) * this.mesh_.ghostCellsm1_ny_bi_[i];

            this.ghostCellsm1_rho_mirror_[i] = rho_mirror;
            this.ghostCellsm1_phi_mirror_[i] = phi_mirror;
            this.ghostCellsm1_u_mirror_[i] = u_mirror;
            this.ghostCellsm1_v_mirror_[i] = v_mirror;
        }

        forall i in 0..<this.ghostCellm1WallIndicesWithGhost_.size {
            const ghostCell = this.ghostCellm1WallIndicesWithGhost_[i];
            this.ghostCellsm1_rho_bi_[i] = 0.5 * (this.rhorho_[ghostCell] + this.ghostCellsm1_rho_mirror_[i]);
            this.ghostCellsm1_phi_bi_[i] = 0.5 * (this.phiphi_[ghostCell] + this.ghostCellsm1_phi_mirror_[i]);
            this.ghostCellsm1_u_bi_[i] = 0.5 * (this.uu_[ghostCell] + this.ghostCellsm1_u_mirror_[i]);
            this.ghostCellsm1_v_bi_[i] = 0.5 * (this.vv_[ghostCell] + this.ghostCellsm1_v_mirror_[i]);
            this.ghostCellsm1_cp_bi_[i] = 2*(this.ghostCellsm1_rho_bi_[i]**this.inputs_.GAMMA_ - 1.0) / (this.inputs_.GAMMA_ * this.inputs_.MACH_**2);
        }


        // Farfield
        if this.inputs_.MACH_ >= 1.0 {
            // Supersonic farfield: TO BE IMPLEMENTED
        }

        else {
            forall face in this.outflowFacesI_ {
                const (i, j) = this.mesh_.IfacesIJ_[face];
                const leftCellWithGhost = meshIndex2FVMindex(i-1, j);
                const rightCellWithGhost = leftCellWithGhost + 1;

                // this.uu_[rightCellWithGhost] = 2*(this.inputs_.U_INF_-this.Rc1_[leftCellWithGhost] / (this.inputs_.C_INF_ - this.inputs_.U_INF_)) - this.uu_[leftCellWithGhost];
                // this.vv_[rightCellWithGhost] = 2*(this.inputs_.V_INF_-this.Rc1_[leftCellWithGhost] / (this.inputs_.C_INF_ - this.inputs_.V_INF_)) - this.vv_[leftCellWithGhost];

                this.uu_[rightCellWithGhost] = this.inputs_.U_INF_;
                this.vv_[rightCellWithGhost] = this.inputs_.V_INF_;
                this.rhorho_[rightCellWithGhost] = this.inputs_.RHO_INF_;
                this.phiphi_[rightCellWithGhost] = this.inputs_.U_INF_ * this.xCellsWithGhosts_[rightCellWithGhost] + this.inputs_.V_INF_ * this.yCellsWithGhosts_[rightCellWithGhost];
            }

            forall face in this.inflowFacesI_ {
                const (i, j) = this.mesh_.IfacesIJ_[face];
                const rightCellWithGhost = meshIndex2FVMindex(i, j);
                const leftCellWithGhost = rightCellWithGhost - 1;

                // this.uu_[leftCellWithGhost] = 2*(this.inputs_.U_INF_+this.Rc1_[rightCellWithGhost] / (this.inputs_.C_INF_ + this.inputs_.U_INF_)) - this.uu_[rightCellWithGhost];
                // this.vv_[leftCellWithGhost] = 2*(this.inputs_.V_INF_+this.Rc1_[rightCellWithGhost] / (this.inputs_.C_INF_ + this.inputs_.V_INF_)) - this.vv_[rightCellWithGhost];

                this.uu_[leftCellWithGhost] = this.inputs_.U_INF_;
                this.vv_[leftCellWithGhost] = this.inputs_.V_INF_;
                this.rhorho_[leftCellWithGhost] = this.inputs_.RHO_INF_;
                this.phiphi_[leftCellWithGhost] = this.inputs_.U_INF_ * this.xCellsWithGhosts_[leftCellWithGhost] + this.inputs_.V_INF_ * this.yCellsWithGhosts_[leftCellWithGhost];
            }

            forall i in 0..<this.mesh_.niCell_ { // bottom boundary
                const face = i;
                const topCellWithGhost = meshIndex2FVMindex(i, 0);
                const bottomCellWithGhost = topCellWithGhost - this.niCellWithGhosts_;

                // this.uu_[bottomCellWithGhost] = 2*(this.inputs_.U_INF_-this.Rc1_[topCellWithGhost] / (this.inputs_.C_INF_ - this.inputs_.U_INF_)) - this.uu_[topCellWithGhost];
                // this.vv_[bottomCellWithGhost] = 2*(this.inputs_.V_INF_-this.Rc1_[topCellWithGhost] / (this.inputs_.C_INF_ - this.inputs_.V_INF_)) - this.vv_[topCellWithGhost];

                this.uu_[bottomCellWithGhost] = this.inputs_.U_INF_;
                this.vv_[bottomCellWithGhost] = this.inputs_.V_INF_;
                this.rhorho_[bottomCellWithGhost] = this.inputs_.RHO_INF_;
                this.phiphi_[bottomCellWithGhost] = this.inputs_.U_INF_ * this.xCellsWithGhosts_[bottomCellWithGhost] + this.inputs_.V_INF_ * this.yCellsWithGhosts_[bottomCellWithGhost];
                
            }

            forall i in 0..<this.mesh_.niCell_ { // top boundary
                const face = i + this.mesh_.niCell_ * this.mesh_.njCell_;
                const bottomCellWithGhost = meshIndex2FVMindex(i, this.mesh_.njCell_ - 1);
                const topCellWithGhost = bottomCellWithGhost + this.niCellWithGhosts_;

                // this.uu_[topCellWithGhost] = 2*(this.inputs_.U_INF_+this.Rc1_[bottomCellWithGhost] / (this.inputs_.C_INF_ + this.inputs_.U_INF_)) - this.uu_[bottomCellWithGhost];
                // this.vv_[topCellWithGhost] = 2*(this.inputs_.V_INF_+this.Rc1_[bottomCellWithGhost] / (this.inputs_.C_INF_ + this.inputs_.V_INF_)) - this.vv_[bottomCellWithGhost];

                this.uu_[topCellWithGhost] = this.inputs_.U_INF_;
                this.vv_[topCellWithGhost] = this.inputs_.V_INF_;
                this.rhorho_[topCellWithGhost] = this.inputs_.RHO_INF_;
                this.phiphi_[topCellWithGhost] = this.inputs_.U_INF_ * this.xCellsWithGhosts_[topCellWithGhost] + this.inputs_.V_INF_ * this.yCellsWithGhosts_[topCellWithGhost];
            }

            // corners
            // bottom left
            {
                const cellWithGhost = 1 + 1 * this.niCellWithGhosts_;
                this.uu_[cellWithGhost] = this.inputs_.U_INF_;
                this.vv_[cellWithGhost] = this.inputs_.V_INF_;
                this.rhorho_[cellWithGhost] = this.inputs_.RHO_INF_;
                this.phiphi_[cellWithGhost] = this.inputs_.U_INF_ * this.xCellsWithGhosts_[cellWithGhost] + this.inputs_.V_INF_ * this.yCellsWithGhosts_[cellWithGhost];
            }
            // bottom right
            {
                const cellWithGhost = (this.niCellWithGhosts_ - 2) + 1 * this.niCellWithGhosts_;
                this.uu_[cellWithGhost] = this.inputs_.U_INF_;
                this.vv_[cellWithGhost] = this.inputs_.V_INF_;
                this.rhorho_[cellWithGhost] = this.inputs_.RHO_INF_;
                this.phiphi_[cellWithGhost] = this.inputs_.U_INF_ * this.xCellsWithGhosts_[cellWithGhost] + this.inputs_.V_INF_ * this.yCellsWithGhosts_[cellWithGhost];
            }
            // top left
            {
                const cellWithGhost = 1 + (this.njCellWithGhosts_ - 2) * this.niCellWithGhosts_;
                this.uu_[cellWithGhost] = this.inputs_.U_INF_;
                this.vv_[cellWithGhost] = this.inputs_.V_INF_;
                this.rhorho_[cellWithGhost] = this.inputs_.RHO_INF_;
                this.phiphi_[cellWithGhost] = this.inputs_.U_INF_ * this.xCellsWithGhosts_[cellWithGhost] + this.inputs_.V_INF_ * this.yCellsWithGhosts_[cellWithGhost];
            }
            // top right
            {
                const cellWithGhost = (this.niCellWithGhosts_ - 2) + (this.njCellWithGhosts_ - 2) * this.niCellWithGhosts_;
                this.uu_[cellWithGhost] = this.inputs_.U_INF_;
                this.vv_[cellWithGhost] = this.inputs_.V_INF_;
                this.rhorho_[cellWithGhost] = this.inputs_.RHO_INF_;
                this.phiphi_[cellWithGhost] = this.inputs_.U_INF_ * this.xCellsWithGhosts_[cellWithGhost] + this.inputs_.V_INF_ * this.yCellsWithGhosts_[cellWithGhost];
            }
            
        }
    }

    proc updatePrimitiveVariablesFromConserved() {
        forall cell in 0..<this.nCellWithGhosts_{
            if this.cellTypesWithGhosts_[cell] != 1 {
                continue;
            }
            const leftCell = cell - 1;
            const rightCell = cell + 1;
            const bottomCell = cell - this.niCellWithGhosts_;
            const topCell = cell + this.niCellWithGhosts_;

            const m_i = this.mi_[cell];
            const n_i = this.ni_[cell];
            const s_i = this.si_[cell];

            const m_j = this.mj_[cell];
            const n_j = this.nj_[cell];
            const s_j = this.sj_[cell];

            this.uu_[cell] = m_i * this.phiphi_[leftCell] + s_i * this.phiphi_[cell] + n_i * this.phiphi_[rightCell];
            this.vv_[cell] = m_j * this.phiphi_[bottomCell] + s_j * this.phiphi_[cell] + n_j * this.phiphi_[topCell];
            this.rhorho_[cell] = (1 + (this.inputs_.GAMMA_ - 1) / 2 * this.inputs_.MACH_**2 * (1 - (this.uu_[cell]**2 + this.vv_[cell]**2)))**this.one_over_gamma_minus_one_;

            // this.uu_[cell] = (this.phiphi_[rightCell] - this.phiphi_[leftCell]) / (this.xCellsWithGhosts_[rightCell] - this.xCellsWithGhosts_[leftCell]);
            // this.vv_[cell] = (this.phiphi_[topCell] - this.phiphi_[bottomCell]) / (this.yCellsWithGhosts_[topCell] - this.yCellsWithGhosts_[bottomCell]);

            // if cell == 81 {
            //     writeln("Cell ", cell, ": u = ", this.uu_[cell], ", v = ", this.vv_[cell], ", Rc0 = ", this.Rc0_[cell], ", Rc1 = ", this.Rc1_[cell]);
            //     writeln("  phi_left = ", this.phiphi_[leftCell], ", phi_center = ", this.phiphi_[cell], ", phi_right = ", this.phiphi_[rightCell]);
            //     writeln("  phi_bottom = ", this.phiphi_[bottomCell], ", phi_center = ", this.phiphi_[cell], ", phi_top = ", this.phiphi_[topCell]);
            //     writeln("  m_i = ", m_i, ", s_i = ", s_i, ", n_i = ", n_i, " sum = ", (m_i + s_i + n_i));
            //     writeln("  m_j = ", m_j, ", s_j = ", s_j, ", n_j = ", n_j, " sum = ", (m_j + s_j + n_j) );
            // }
        }
    }

    proc compute_convective_fluxes() {
        forall face in 0..<this.mesh_.nFace_ {
            const (i, j) = this.mesh_.IfacesIJ_[face];
            const leftCell = meshIndex2FVMindex(i-1, j);
            const rightCell = meshIndex2FVMindex(i, j);

            if ((this.cellTypesWithGhosts_[leftCell] != 1 && this.cellTypesWithGhosts_[rightCell] != 1) 
            && (this.cellTypesWithGhosts_[leftCell] != 0 && this.cellTypesWithGhosts_[rightCell] != 0)) {
                continue;
            }

            // const avg_rho = 0.5 * (this.rhorho_[leftCell] + this.rhorho_[rightCell]);
            // const avg_u = 0.5 * (this.uu_[leftCell] + this.uu_[rightCell]);
            // const avg_v = 0.5 * (this.vv_[leftCell] + this.vv_[rightCell]);

            const avg_u = (this.phiphi_[rightCell] - this.phiphi_[leftCell]) / (this.xCellsWithGhosts_[rightCell] - this.xCellsWithGhosts_[leftCell]);
            const v_right = (this.phiphi_[rightCell + this.niCellWithGhosts_] - this.phiphi_[rightCell - this.niCellWithGhosts_]) / (this.yCellsWithGhosts_[rightCell + this.niCellWithGhosts_] - this.yCellsWithGhosts_[rightCell - this.niCellWithGhosts_]);
            const v_left = (this.phiphi_[leftCell + this.niCellWithGhosts_] - this.phiphi_[leftCell - this.niCellWithGhosts_]) / (this.yCellsWithGhosts_[leftCell + this.niCellWithGhosts_] - this.yCellsWithGhosts_[leftCell - this.niCellWithGhosts_]);
            const avg_v = 0.5 * (v_left + v_right);
            // const avg_rho = (1 + (this.inputs_.GAMMA_ - 1) / 2 * this.inputs_.MACH_**2 * (1 - (avg_u**2 + avg_v**2)))**this.one_over_gamma_minus_one_;

            const avg_rho = (1 + (this.inputs_.GAMMA_ - 1) / 2 * this.inputs_.MACH_**2 * (1 - (avg_u**2 + avg_v**2) - this.phi_t_[leftCell] - this.phi_t_[rightCell]))**this.one_over_gamma_minus_one_;

            this.F0I_[face] = avg_rho * avg_u * this.mesh_.IfaceAreas_[face];
            this.dF0Idrho_[face] = avg_u * this.mesh_.IfaceAreas_[face];
            this.avg_rho_I_[face] = avg_rho;
        }

        forall face in 0..<this.mesh_.nFace_ {
            const (i, j) = this.mesh_.JfacesIJ_[face];
            const bottomCell = meshIndex2FVMindex(i, j-1);
            const topCell = meshIndex2FVMindex(i, j);

            if ((this.cellTypesWithGhosts_[bottomCell] != 1 && this.cellTypesWithGhosts_[topCell] != 1) 
            && (this.cellTypesWithGhosts_[bottomCell] != 0 && this.cellTypesWithGhosts_[topCell] != 0)) {
                continue;
            }

            // const avg_rho = 0.5 * (this.rhorho_[bottomCell] + this.rhorho_[topCell]);
            // const avg_u = 0.5 * (this.uu_[bottomCell] + this.uu_[topCell]);
            // const avg_v = 0.5 * (this.vv_[bottomCell] + this.vv_[topCell]);

            const avg_v = (this.phiphi_[topCell] - this.phiphi_[bottomCell]) / (this.yCellsWithGhosts_[topCell] - this.yCellsWithGhosts_[bottomCell]);
            const u_top = (this.phiphi_[topCell + 1] - this.phiphi_[topCell - 1]) / (this.xCellsWithGhosts_[topCell + 1] - this.xCellsWithGhosts_[topCell - 1]);
            const u_bottom = (this.phiphi_[bottomCell + 1] - this.phiphi_[bottomCell - 1]) / (this.xCellsWithGhosts_[bottomCell + 1] - this.xCellsWithGhosts_[bottomCell - 1]);
            const avg_u = 0.5 * (u_bottom + u_top);
            // const avg_rho = (1 + (this.inputs_.GAMMA_ - 1) / 2 * this.inputs_.MACH_**2 * (1 - (avg_u**2 + avg_v**2)))**this.one_over_gamma_minus_one_;

            const avg_rho = (1 + (this.inputs_.GAMMA_ - 1) / 2 * this.inputs_.MACH_**2 * (1 - (avg_u**2 + avg_v**2) - this.phi_t_[bottomCell] - this.phi_t_[topCell]))**this.one_over_gamma_minus_one_;

            this.F0J_[face] = avg_rho * avg_v * this.mesh_.JfaceAreas_[face];
            this.dF0Jdrho_[face] = avg_v * this.mesh_.JfaceAreas_[face];
            this.avg_rho_J_[face] = avg_rho;
        }
    }

    proc compute_lambdas() {
        forall j in 0..<this.mesh_.njCell_ {
            for i in 0..<this.mesh_.niCell_ {
                const idx = this.mesh_.iiCell(i, j);
                if this.mesh_.cellTypes_[idx] != 1 && this.mesh_.cellTypes_[idx] != 0 && this.mesh_.cellTypes_[idx] != -2 {
                    continue;
                }
                const idxWithGhost = meshIndex2FVMindex(i, j);
                const u = this.uu_[idxWithGhost];
                const v = this.vv_[idxWithGhost];
                const rho = this.rhorho_[idxWithGhost];
                const a = sqrt(rho**(this.inputs_.GAMMA_ - 1) / this.inputs_.MACH_**2);

                this.LambdaI_[idxWithGhost] = (abs(u) + a) * this.mesh_.avgFaceAreaI_[idx];
                this.LambdaJ_[idxWithGhost] = (abs(v) + a) * this.mesh_.avgFaceAreaJ_[idx];
            }
        }
    }

    proc epsilon(p_Im1: real(64), p_I: real(64), p_Ip1: real(64), p_Ip2: real(64)) {
        const gamma_I = abs(p_Ip1 - 2.0 * p_I + p_Im1) / (p_Ip1 + 2.0 * p_I + p_Im1);
        const gamma_Ip1 = abs(p_Ip2 - 2.0 * p_Ip1 + p_I) / (p_Ip2 + 2.0 * p_Ip1 + p_I);
        const eps2 = this.inputs_.K2_ * max(gamma_I, gamma_Ip1);
        const eps4 = max(0.0, this.inputs_.K4_ - eps2);
        return (eps2, eps4);
    }

    proc compute_diffusive_fluxes() {
        forall j in 0..<this.mesh_.njCell_ {
            for i in 0..<this.mesh_.niCell_ {
                const idxWithGhost = meshIndex2FVMindex(i, j);
                this.pp_[idxWithGhost] = this.rhorho_[idxWithGhost]**this.inputs_.GAMMA_ / (this.inputs_.GAMMA_ * this.inputs_.MACH_**2);
            }
        }

        forall face in 0..<this.mesh_.nFace_ {
            const (i, j) = this.mesh_.IfacesIJ_[face];
            const leftCell = meshIndex2FVMindex(i-1, j);
            const rightCell = meshIndex2FVMindex(i, j);

            if (this.cellTypesWithGhosts_[leftCell] != 1 || this.cellTypesWithGhosts_[rightCell] != 1) {
                continue;
            }

            const leftCell_m1 = leftCell - 1;
            const rightCell_p1 = rightCell + 1;

            const (eps2, eps4) = this.epsilon(this.pp_[leftCell_m1], this.pp_[leftCell], this.pp_[rightCell], this.pp_[rightCell_p1]);
            const LambdaS = 0.5*(this.LambdaI_[leftCell] + this.LambdaI_[rightCell]) + 0.5*(this.LambdaJ_[leftCell] + this.LambdaJ_[rightCell]);

            this.D0I_[face] = LambdaS * (eps2 * (this.rhorho_[rightCell] - this.rhorho_[leftCell]) - eps4 * (this.rhorho_[rightCell_p1] - 3.0 * this.rhorho_[rightCell] + 3.0 * this.rhorho_[leftCell] - this.rhorho_[leftCell_m1]) );
            
        }

        forall face in 0..<this.mesh_.nFace_ {
            const (i, j) = this.mesh_.JfacesIJ_[face];
            const bottomCell = meshIndex2FVMindex(i, j-1);
            const topCell = meshIndex2FVMindex(i, j);

            if (this.cellTypesWithGhosts_[bottomCell] != 1 || this.cellTypesWithGhosts_[topCell] != 1) {
                continue;
            }

            const bottomCell_m1 = bottomCell - this.niCellWithGhosts_;
            const topCell_p1 = topCell + this.niCellWithGhosts_;

            const (eps2, eps4) = this.epsilon(this.pp_[bottomCell_m1], this.pp_[bottomCell], this.pp_[topCell], this.pp_[topCell_p1]);
            const LambdaS = 0.5*(this.LambdaI_[bottomCell] + this.LambdaI_[topCell]) + 0.5*(this.LambdaJ_[bottomCell] + this.LambdaJ_[topCell]);

            this.D0J_[face] = LambdaS * (eps2 * (this.rhorho_[topCell] - this.rhorho_[bottomCell]) - eps4 * (this.rhorho_[topCell_p1] - 3.0 * this.rhorho_[topCell] + 3.0 * this.rhorho_[bottomCell] - this.rhorho_[bottomCell_m1]) );
            
        }
    }

    proc compute_convective_residuals() {
        this.Rc0_ = 0.0;
        this.Rc1_ = 0.0;
        
        forall j in 0..<this.mesh_.njCell_ {
            for i in 0..<this.mesh_.niCell_ {
                const idxWithGhost = meshIndex2FVMindex(i, j);
                if this.cellTypesWithGhosts_[idxWithGhost] != 1 {
                    continue;
                }
                const leftFace = i + j*this.mesh_.niNode_;
                const rightFace = leftFace + 1;
                const bottomFace = i + j*this.mesh_.niCell_;
                const topFace = bottomFace + this.mesh_.niCell_;

                this.Rc0_[idxWithGhost] = -F0I_[leftFace] + F0I_[rightFace] - F0J_[bottomFace] + F0J_[topFace];
                this.Rc1_[idxWithGhost] = this.cellVolumesWithGhosts_[idxWithGhost]*(0.5*(this.uu_[idxWithGhost]**2 + this.vv_[idxWithGhost]**2) + 
                                        ((this.rhorho_[idxWithGhost]**(this.inputs_.GAMMA_ - 1.0) - 1.0) /
                                        (this.inputs_.MACH_**2 * (this.inputs_.GAMMA_ - 1.0)) - 0.5));
                // this.Rc1_[idxWithGhost] = (0.5*(this.uu_[idxWithGhost]**2 + this.vv_[idxWithGhost]**2) + 
                //                         ((this.rhorho_[idxWithGhost]**(this.inputs_.GAMMA_ - 1.0) - 1.0) /
                //                         (this.inputs_.MACH_**2 * (this.inputs_.GAMMA_ - 1.0)) - 0.5));

                // const idx = this.mesh_.iiCell(i, j);
                // if idx == 0 {
                //     writeln("Cell ", idxWithGhost, ": Rc0 = ", this.Rc0_[idxWithGhost], ", Rc1 = ", this.Rc1_[idxWithGhost]);
                //     writeln("  F0I_left = ", F0I_[leftFace], ", F0I_right = ", F0I_[rightFace]);
                //     writeln("  F0J_bottom = ", F0J_[bottomFace], ", F0J_top = ", F0J_[topFace]);
                // }                                     
            }
        }
    }

    proc compute_diffusive_residuals() {
        this.Rd0_ = 0.0;
        this.Rd1_ = 0.0;
        forall j in 0..<this.mesh_.njCell_ {
            for i in 0..<this.mesh_.niCell_ {
                const idxWithGhost = meshIndex2FVMindex(i, j);
                if this.cellTypesWithGhosts_[idxWithGhost] != 1 {
                    continue;
                }
                const leftFace = i + j*this.mesh_.niNode_;
                const rightFace = leftFace + 1;
                const bottomFace = i + j*this.mesh_.niCell_;
                const topFace = bottomFace + this.mesh_.niCell_;

                this.Rd0_[idxWithGhost] = -D0I_[leftFace] + D0I_[rightFace] - D0J_[bottomFace] + D0J_[topFace];
                this.Rd1_[idxWithGhost] = -D1I_[leftFace] + D1I_[rightFace] - D1J_[bottomFace] + D1J_[topFace];
            }
        }
    }

    proc run() {
        this.updateGhostCells();
        this.updatePrimitiveVariablesFromConserved();
        this.compute_lambdas();
        this.compute_convective_fluxes();
        this.compute_convective_residuals();
        // this.compute_diffusive_fluxes();
        // this.compute_diffusive_residuals();
    }

    proc compute_aerodynamics_coefficients() {
        var Fx = 0.0;
        var Fy = 0.0;
        var M = 0.0;

        var avg_nx : [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var avg_ny : [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var avg_cp : [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var avg_x : [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var avg_y : [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var avg_area : [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);

        forall i in 0..<this.ghostCellWallIndicesWithGhost_.size {
            const ip1 = (i + 1) % this.ghostCellWallIndicesWithGhost_.size;
            avg_nx[i] = -0.5 * (this.mesh_.ghostCells_nx_bi_[i] + this.mesh_.ghostCells_nx_bi_[ip1]);
            avg_ny[i] = -0.5 * (this.mesh_.ghostCells_ny_bi_[i] + this.mesh_.ghostCells_ny_bi_[ip1]);
            avg_cp[i] = 0.5 * (this.ghostCells_cp_bi_[i] + this.ghostCells_cp_bi_[ip1]);
            avg_x[i] = 0.5 * (this.mesh_.ghostCells_x_bi_[i] + this.mesh_.ghostCells_x_bi_[ip1]);
            avg_y[i] = 0.5 * (this.mesh_.ghostCells_y_bi_[i] + this.mesh_.ghostCells_y_bi_[ip1]);
            avg_area[i] = sqrt( (this.mesh_.ghostCells_x_bi_[ip1] - this.mesh_.ghostCells_x_bi_[i])**2 + (this.mesh_.ghostCells_y_bi_[ip1] - this.mesh_.ghostCells_y_bi_[i])**2 );
        }

        forall i in 0..<this.ghostCellWallIndicesWithGhost_.size with (+reduce Fx, +reduce Fy, +reduce M) {
            Fx += avg_cp[i] * avg_nx[i] * avg_area[i];
            Fy += avg_cp[i] * avg_ny[i] * avg_area[i];
            M += avg_cp[i] * ((this.inputs_.X_REF_ - avg_x[i]) * avg_ny[i] + (avg_y[i] - this.inputs_.Y_REF_) * avg_nx[i]) * avg_area[i];
        }

        const Cl = Fy*cos(this.inputs_.ALPHA_ * pi / 180.0) - Fx*sin(this.inputs_.ALPHA_ * pi / 180.0);
        const Cd = Fx*cos(this.inputs_.ALPHA_ * pi / 180.0) + Fy*sin(this.inputs_.ALPHA_ * pi / 180.0);

        return (Cl, Cd, M);
    }

    proc meshIndex2FVMindex(i: int, j: int): int {
        return i+2 + (j+2) * this.niCellWithGhosts_;
    }

    proc writeSolution2CGNS(time: list(real(64)),
                                 iterations: list(int),
                                 res0: list(real(64)),
                                 res1: list(real(64)),
                                 cls: list(real(64)),
                                 cds: list(real(64)),
                                 cms: list(real(64))) {
        const dom = this.mesh_.cell_dom_;
        var rho: [dom] real(64);
        var u: [dom] real(64);
        var v: [dom] real(64);
        var w: [dom] real(64);
        var p: [dom] real(64);
        var Rc0: [dom] real(64);
        var Rc1: [dom] real(64);
        var Rd0: [dom] real(64);
        var Rd1: [dom] real(64);
        var R0: [dom] real(64);
        var R1: [dom] real(64);
        var mach: [dom] real(64);
        var phi: [dom] real(64);

        forall j in 0..<this.mesh_.njCell_ {
            for i in 0..<this.mesh_.niCell_ {
                const idxWithGhost = meshIndex2FVMindex(i, j);
                this.pp_[idxWithGhost] = this.rhorho_[idxWithGhost]**this.inputs_.GAMMA_ / (this.inputs_.GAMMA_ * this.inputs_.MACH_**2);
            }
        }

        forall j in 0..<this.mesh_.njCell_ {
            for i in 0..<this.mesh_.niCell_ {
                const idx = this.mesh_.iiCell(i, j);
                const idxWithGhost = meshIndex2FVMindex(i, j);
                rho[idx] = this.rhorho_[idxWithGhost];
                u[idx] = this.uu_[idxWithGhost];
                v[idx] = this.vv_[idxWithGhost];
                p[idx] = this.pp_[idxWithGhost];
                
                Rc0[idx] = this.Rc0_[idxWithGhost];
                Rc1[idx] = this.Rc1_[idxWithGhost];
                Rd0[idx] = this.Rd0_[idxWithGhost];
                Rd1[idx] = this.Rd1_[idxWithGhost];
                R0[idx] = this.R0_[idxWithGhost];
                R1[idx] = this.R1_[idxWithGhost];
                phi[idx] = this.phiphi_[idxWithGhost];
                mach[idx] = sqrt(u[idx]**2 + v[idx]**2) * this.inputs_.MACH_ * rho[idx]**((1.0 - this.inputs_.GAMMA_)/2.0);
            }
        }

        var fields = new map(string, [dom] real(64));
        fields["Density"] = rho;
        fields["VelocityX"] = u;
        fields["VelocityY"] = v;
        fields["VelocityZ"] = w;
        fields["Pressure"] = p;
        fields["Rc0"] = Rc0;
        fields["Rc1"] = Rc1;
        fields["Rd0"] = Rd0;
        fields["Rd1"] = Rd1;
        fields["R0"] = R0;
        fields["R1"] = R1;
        fields["Mach"] = mach;
        fields["Phi"] = phi;



        var writer = new cgnsFlowWriter_c(this.inputs_.OUTPUT_FILENAME_);
        writer.writeToCGNS(this.mesh_, dom, fields);

        writer.writeConvergenceHistory(time, iterations, res0, res1, cls, cds, cms);

        // Wall solution
        var TEkls_ = this.mesh_.TEkls_!.borrow();
        for (j, nearestCell) in zip(0..<4, this.mesh_.nearestCellTE_[0..]) { // MAYBE GENERALIZE LATER
            TEkls_.rho_fieldValues_[j] = this.rhorho_[nearestCell];
            TEkls_.u_fieldValues_[j] = this.uu_[nearestCell];
            TEkls_.v_fieldValues_[j] = this.vv_[nearestCell];
            TEkls_.p_fieldValues_[j] = this.pp_[nearestCell];
        }
        const rho_TE = TEkls_.interpolate(TEkls_.rho_fieldValues_);
        const u_TE = TEkls_.interpolate(TEkls_.u_fieldValues_);
        const v_TE = TEkls_.interpolate(TEkls_.v_fieldValues_);
        const p_TE = TEkls_.interpolate(TEkls_.p_fieldValues_);
        const mach_TE = sqrt(u_TE**2 + v_TE**2) * this.inputs_.MACH_ * rho_TE**((1.0 - this.inputs_.GAMMA_)/2.0);
        const cp_TE = (p_TE - this.inputs_.P_INF_) / (0.5 * this.inputs_.RHO_INF_ * (this.inputs_.U_INF_**2 + this.inputs_.V_INF_**2));

        var LEkls_ = this.mesh_.LEkls_!.borrow();
        for (j, nearestCell) in zip(0..<4, this.mesh_.nearestCellLE_[0..]) { // MAYBE GENERALIZE LATER
            LEkls_.rho_fieldValues_[j] = this.rhorho_[nearestCell];
            LEkls_.u_fieldValues_[j] = this.uu_[nearestCell];
            LEkls_.v_fieldValues_[j] = this.vv_[nearestCell];
            LEkls_.p_fieldValues_[j] = this.pp_[nearestCell];
        }
        const rho_LE = LEkls_.interpolate(LEkls_.rho_fieldValues_);
        const u_LE = LEkls_.interpolate(LEkls_.u_fieldValues_);
        const v_LE = LEkls_.interpolate(LEkls_.v_fieldValues_);
        const p_LE = LEkls_.interpolate(LEkls_.p_fieldValues_);
        const mach_LE = sqrt(u_LE**2 + v_LE**2) * this.inputs_.MACH_ * rho_LE**((1.0 - this.inputs_.GAMMA_)/2.0);
        const cp_LE = (p_LE - this.inputs_.P_INF_) / (0.5 * this.inputs_.RHO_INF_ * (this.inputs_.U_INF_**2 + this.inputs_.V_INF_**2));

        var wall_fields_LE_TE = new map(string, [0..<2] real(64));
        wall_fields_LE_TE["Density"] = [rho_LE, rho_TE];
        wall_fields_LE_TE["VelocityX"] = [u_LE, u_TE];
        wall_fields_LE_TE["VelocityY"] = [v_LE, v_TE];
        wall_fields_LE_TE["Pressure"] = [p_LE, p_TE];
        wall_fields_LE_TE["Mach"] = [mach_LE, mach_TE];
        wall_fields_LE_TE["Cp"] = [cp_LE, cp_TE];

        writer.writeLEandTEtoCGNS(this.mesh_, {0..<2},wall_fields_LE_TE);


        var wall_mach: [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var ghostCx: [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var ghostCy: [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var ghostuu: [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var ghostvv: [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var ghostpp: [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var ghostrho: [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        forall i in 0..<this.ghostCellWallIndicesWithGhost_.size {
            const idx = this.mesh_.ghostCellIndices_[i];
            const idxWithGhost = this.ghostCellWallIndicesWithGhost_[i];
            const rho_w = this.ghostCells_rho_bi_[i];
            const u_w = this.ghostCells_u_bi_[i];
            const v_w = this.ghostCells_v_bi_[i];
            wall_mach[i] = sqrt(u_w**2 + v_w**2) * this.inputs_.MACH_ * rho_w**((1.0 - this.inputs_.GAMMA_)/2.0);
            ghostCx[i] = this.mesh_.xCells_[idx];
            ghostCy[i] = this.mesh_.yCells_[idx];
            ghostuu[i] = this.uu_[idxWithGhost];
            ghostvv[i] = this.vv_[idxWithGhost];
            ghostpp[i] = this.pp_[idxWithGhost];
            ghostrho[i] = this.rhorho_[idxWithGhost];
        }

        var wall_fields = new map(string, [0..<this.ghostCellWallIndicesWithGhost_.size] real(64));
        wall_fields["x"] = this.mesh_.ghostCells_x_bi_;
        wall_fields["y"] = this.mesh_.ghostCells_y_bi_;
        wall_fields["x_mirror"] = this.mesh_.ghostCells_x_mirror_;
        wall_fields["y_mirror"] = this.mesh_.ghostCells_y_mirror_;
        wall_fields["Density"] = this.ghostCells_rho_bi_;
        wall_fields["VelocityX"] = this.ghostCells_u_bi_;
        wall_fields["VelocityY"] = this.ghostCells_v_bi_;
        wall_fields["VelocityZ"] = this.ghostCells_w_bi_;
        wall_fields["Cp"] = this.ghostCells_cp_bi_;
        wall_fields["nx"] = this.mesh_.ghostCells_nx_bi_;
        wall_fields["ny"] = this.mesh_.ghostCells_ny_bi_;
        wall_fields["curvature"] = this.mesh_.ghostCells_curvature_bi_;
        wall_fields["Mach"] = wall_mach;
        wall_fields["x_ghost"] = ghostCx;
        wall_fields["y_ghost"] = ghostCy;
        wall_fields["VelocityX_ghost"] = ghostuu;
        wall_fields["VelocityY_ghost"] = ghostvv;
        wall_fields["Pressure_ghost"] = ghostpp;
        wall_fields["Density_ghost"] = ghostrho;

        writer.writeWallToCGNS(this.mesh_, {0..<this.ghostCellWallIndicesWithGhost_.size},wall_fields);

        // wall2 solution
        var wallm1_mach : [0..<this.ghostCellm1WallIndicesWithGhost_.size] real(64);
        var ghostm1Cx : [0..<this.ghostCellm1WallIndicesWithGhost_.size] real(64);
        var ghostm1Cy : [0..<this.ghostCellm1WallIndicesWithGhost_.size] real(64);
        var ghostm1uu : [0..<this.ghostCellm1WallIndicesWithGhost_.size] real(64);
        var ghostm1vv : [0..<this.ghostCellm1WallIndicesWithGhost_.size] real(64);
        var ghostm1pp : [0..<this.ghostCellm1WallIndicesWithGhost_.size] real(64);
        var ghostm1rho : [0..<this.ghostCellm1WallIndicesWithGhost_.size] real(64);
        forall i in 0..<this.ghostCellm1WallIndicesWithGhost_.size {
            const idx = this.mesh_.ghostCellm1Indices_[i];
            const idxWithGhost = this.ghostCellm1WallIndicesWithGhost_[i];
            const rho_w = this.ghostCellsm1_rho_bi_[i];
            const u_w = this.ghostCellsm1_u_bi_[i];
            const v_w = this.ghostCellsm1_v_bi_[i];
            wallm1_mach[i] = sqrt(u_w**2 + v_w**2) * this.inputs_.MACH_ * rho_w**((1.0 - this.inputs_.GAMMA_)/2.0);
            ghostm1Cx[i] = this.mesh_.xCells_[idx];
            ghostm1Cy[i] = this.mesh_.yCells_[idx];
            ghostm1uu[i] = this.uu_[idxWithGhost];
            ghostm1vv[i] = this.vv_[idxWithGhost];
            ghostm1pp[i] = this.pp_[idxWithGhost];
            ghostm1rho[i] = this.rhorho_[idxWithGhost];
        }

        var wallm1_fields = new map(string, [0..<this.ghostCellm1WallIndicesWithGhost_.size] real(64));
        wallm1_fields["x"] = this.mesh_.ghostCellsm1_x_bi_;
        wallm1_fields["y"] = this.mesh_.ghostCellsm1_y_bi_;
        wallm1_fields["x_mirror"] = this.mesh_.ghostCellsm1_x_mirror_;
        wallm1_fields["y_mirror"] = this.mesh_.ghostCellsm1_y_mirror_;
        wallm1_fields["Density"] = this.ghostCellsm1_rho_bi_;
        wallm1_fields["VelocityX"] = this.ghostCellsm1_u_bi_;
        wallm1_fields["VelocityY"] = this.ghostCellsm1_v_bi_;
        wallm1_fields["VelocityZ"] = this.ghostCellsm1_w_bi_;
        wallm1_fields["Cp"] = this.ghostCellsm1_cp_bi_;
        wallm1_fields["nx"] = this.mesh_.ghostCellsm1_nx_bi_;
        wallm1_fields["ny"] = this.mesh_.ghostCellsm1_ny_bi_;
        wallm1_fields["Mach"] = wallm1_mach;
        wallm1_fields["x_ghost"] = ghostm1Cx;
        wallm1_fields["y_ghost"] = ghostm1Cy;
        wallm1_fields["VelocityX_ghost"] = ghostm1uu;;
        wallm1_fields["VelocityY_ghost"] = ghostm1vv;
        wallm1_fields["Pressure_ghost"] = ghostm1pp;
        wallm1_fields["Density_ghost"] = ghostm1rho;

        writer.writeWall2ToCGNS(this.mesh_, {0..<this.ghostCellm1WallIndicesWithGhost_.size},wallm1_fields);

    }

    proc writeSolution2CGNS() {
        const dom = this.mesh_.cell_dom_;
        var rho: [dom] real(64);
        var u: [dom] real(64);
        var v: [dom] real(64);
        var w: [dom] real(64);
        var p: [dom] real(64);
        var Rc0: [dom] real(64);
        var Rc1: [dom] real(64);
        var Rd0: [dom] real(64);
        var Rd1: [dom] real(64);
        var R0: [dom] real(64);
        var R1: [dom] real(64);
        var mach: [dom] real(64);

        forall j in 0..<this.mesh_.njCell_ {
            for i in 0..<this.mesh_.niCell_ {
                const idxWithGhost = meshIndex2FVMindex(i, j);
                this.pp_[idxWithGhost] = this.rhorho_[idxWithGhost]**this.inputs_.GAMMA_ / (this.inputs_.GAMMA_ * this.inputs_.MACH_**2);
            }
        }

        forall j in 0..<this.mesh_.njCell_ {
            for i in 0..<this.mesh_.niCell_ {
                const idx = this.mesh_.iiCell(i, j);
                const idxWithGhost = meshIndex2FVMindex(i, j);
                rho[idx] = this.rhorho_[idxWithGhost];
                u[idx] = this.uu_[idxWithGhost];
                v[idx] = this.vv_[idxWithGhost];
                p[idx] = this.pp_[idxWithGhost];

                if (this.cellTypesWithGhosts_[idxWithGhost] == 1) {
                    Rc0[idx] = this.Rc0_[idxWithGhost];
                    Rc1[idx] = this.Rc1_[idxWithGhost];
                    Rd0[idx] = this.Rd0_[idxWithGhost];
                    Rd1[idx] = this.Rd1_[idxWithGhost];
                    R0[idx] = this.R0_[idxWithGhost];
                    R1[idx] = this.R1_[idxWithGhost];
                }
                mach[idx] = sqrt(u[idx]**2 + v[idx]**2) * this.inputs_.MACH_ * rho[idx]**((1.0 - this.inputs_.GAMMA_)/2.0);
            }
        }

        var fields = new map(string, [dom] real(64));
        fields["Density"] = rho;
        fields["VelocityX"] = u;
        fields["VelocityY"] = v;
        fields["VelocityZ"] = w;
        fields["Pressure"] = p;
        fields["Rc0"] = Rc0;
        fields["Rc1"] = Rc1;
        fields["Rd0"] = Rd0;
        fields["Rd1"] = Rd1;
        fields["R0"] = R0;
        fields["R1"] = R1;
        fields["Mach"] = mach;



        var writer = new cgnsFlowWriter_c(this.inputs_.OUTPUT_FILENAME_);
        writer.writeToCGNS(this.mesh_, dom, fields);

        // Wall solution
        var TEkls_ = this.mesh_.TEkls_!.borrow();
        for (j, nearestCell) in zip(0..<4, this.mesh_.nearestCellTE_[0..]) { // MAYBE GENERALIZE LATER
            TEkls_.rho_fieldValues_[j] = this.rhorho_[nearestCell];
            TEkls_.u_fieldValues_[j] = this.uu_[nearestCell];
            TEkls_.v_fieldValues_[j] = this.vv_[nearestCell];
            TEkls_.p_fieldValues_[j] = this.pp_[nearestCell];
        }
        const rho_TE = TEkls_.interpolate(TEkls_.rho_fieldValues_);
        const u_TE = TEkls_.interpolate(TEkls_.u_fieldValues_);
        const v_TE = TEkls_.interpolate(TEkls_.v_fieldValues_);
        const p_TE = TEkls_.interpolate(TEkls_.p_fieldValues_);
        const mach_TE = sqrt(u_TE**2 + v_TE**2) * this.inputs_.MACH_ * rho_TE**((1.0 - this.inputs_.GAMMA_)/2.0);
        const cp_TE = (p_TE - this.inputs_.P_INF_) / (0.5 * this.inputs_.RHO_INF_ * (this.inputs_.U_INF_**2 + this.inputs_.V_INF_**2));

        var LEkls_ = this.mesh_.LEkls_!.borrow();
        for (j, nearestCell) in zip(0..<4, this.mesh_.nearestCellLE_[0..]) { // MAYBE GENERALIZE LATER
            LEkls_.rho_fieldValues_[j] = this.rhorho_[nearestCell];
            LEkls_.u_fieldValues_[j] = this.uu_[nearestCell];
            LEkls_.v_fieldValues_[j] = this.vv_[nearestCell];
            LEkls_.p_fieldValues_[j] = this.pp_[nearestCell];
        }
        const rho_LE = LEkls_.interpolate(LEkls_.rho_fieldValues_);
        const u_LE = LEkls_.interpolate(LEkls_.u_fieldValues_);
        const v_LE = LEkls_.interpolate(LEkls_.v_fieldValues_);
        const p_LE = LEkls_.interpolate(LEkls_.p_fieldValues_);
        const mach_LE = sqrt(u_LE**2 + v_LE**2) * this.inputs_.MACH_ * rho_LE**((1.0 - this.inputs_.GAMMA_)/2.0);
        const cp_LE = (p_LE - this.inputs_.P_INF_) / (0.5 * this.inputs_.RHO_INF_ * (this.inputs_.U_INF_**2 + this.inputs_.V_INF_**2));

        var wall_fields_LE_TE = new map(string, [0..<2] real(64));
        wall_fields_LE_TE["Density"] = [rho_LE, rho_TE];
        wall_fields_LE_TE["VelocityX"] = [u_LE, u_TE];
        wall_fields_LE_TE["VelocityY"] = [v_LE, v_TE];
        wall_fields_LE_TE["Pressure"] = [p_LE, p_TE];
        wall_fields_LE_TE["Mach"] = [mach_LE, mach_TE];
        wall_fields_LE_TE["Cp"] = [cp_LE, cp_TE];

        writer.writeLEandTEtoCGNS(this.mesh_, {0..<2},wall_fields_LE_TE);


        var wall_mach: [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var ghostCx: [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var ghostCy: [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var ghostuu: [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var ghostvv: [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var ghostpp: [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        var ghostrho: [0..<this.ghostCellWallIndicesWithGhost_.size] real(64);
        forall i in 0..<this.ghostCellWallIndicesWithGhost_.size {
            const idx = this.mesh_.ghostCellIndices_[i];
            const idxWithGhost = this.ghostCellWallIndicesWithGhost_[i];
            const rho_w = this.ghostCells_rho_bi_[i];
            const u_w = this.ghostCells_u_bi_[i];
            const v_w = this.ghostCells_v_bi_[i];
            wall_mach[i] = sqrt(u_w**2 + v_w**2) * this.inputs_.MACH_ * rho_w**((1.0 - this.inputs_.GAMMA_)/2.0);
            ghostCx[i] = this.mesh_.xCells_[idx];
            ghostCy[i] = this.mesh_.yCells_[idx];
            ghostuu[i] = this.uu_[idxWithGhost];
            ghostvv[i] = this.vv_[idxWithGhost];
            ghostpp[i] = this.pp_[idxWithGhost];
            ghostrho[i] = this.rhorho_[idxWithGhost];
        }

        var wall_fields = new map(string, [0..<this.ghostCellWallIndicesWithGhost_.size] real(64));
        wall_fields["x"] = this.mesh_.ghostCells_x_bi_;
        wall_fields["y"] = this.mesh_.ghostCells_y_bi_;
        wall_fields["x_mirror"] = this.mesh_.ghostCells_x_mirror_;
        wall_fields["y_mirror"] = this.mesh_.ghostCells_y_mirror_;
        wall_fields["Density"] = this.ghostCells_rho_bi_;
        wall_fields["VelocityX"] = this.ghostCells_u_bi_;
        wall_fields["VelocityY"] = this.ghostCells_v_bi_;
        wall_fields["VelocityZ"] = this.ghostCells_w_bi_;
        wall_fields["Cp"] = this.ghostCells_cp_bi_;
        wall_fields["nx"] = this.mesh_.ghostCells_nx_bi_;
        wall_fields["ny"] = this.mesh_.ghostCells_ny_bi_;
        wall_fields["curvature"] = this.mesh_.ghostCells_curvature_bi_;
        wall_fields["Mach"] = wall_mach;
        wall_fields["x_ghost"] = ghostCx;
        wall_fields["y_ghost"] = ghostCy;
        wall_fields["VelocityX_ghost"] = ghostuu;
        wall_fields["VelocityY_ghost"] = ghostvv;
        wall_fields["Pressure_ghost"] = ghostpp;
        wall_fields["Density_ghost"] = ghostrho;

        writer.writeWallToCGNS(this.mesh_, {0..<this.ghostCellWallIndicesWithGhost_.size},wall_fields);

    }

    
}





} // module fullPotentialSpatialDiscretization