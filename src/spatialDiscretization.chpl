module spatialDiscretization 
{
use mesh;
use Time;
use Random;
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

    var Rd0_ : [cell_dom_with_ghosts_] real(64);
    var Rd1_ : [cell_dom_with_ghosts_] real(64);
    var Rd2_ : [cell_dom_with_ghosts_] real(64);
    var Rd3_ : [cell_dom_with_ghosts_] real(64);

    var ghostCellIndices_dom_ : domain(1) = {1..0};
    var ghostCellIndices_ : [ghostCellIndices_dom_] int;

    var cell_dom_ : domain(1) = {1..0};
    var LambdaI_ : [cell_dom_] real(64);
    var LambdaJ_ : [cell_dom_] real(64);

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
        this.cell_dom_ = {0..<(niCell_ * njCell_)};
        this.face_dom_ = {0..<( (niCell_+1)*njCell_ )};
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
        var id = 0;
        for (i, j) in this.mesh_.ghostCellIJ_ {
            const idxWithGhost = i+2 + (j+2) * this.niCellWithGhosts_;
            this.ghostCellIndices_[id] = idxWithGhost;
            id += 1;
        }
    }

    proc updateGhostCells() {
        // Wall
        


        // Farfield
    }

}





}