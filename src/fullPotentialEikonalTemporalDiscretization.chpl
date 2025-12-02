module fullPotentialEikonalTemporalDiscretization 
{
use mesh;
use fullPotentialEikonalSpatialDiscretization;
use linearAlgebra;
import input.inputsConfig;
use Time;
use Math;
use List;

class fullPotentialEikonalTemporalDiscretization {
    var spatialDisc_ : shared fullPotentialEikonalSpatialDiscretization;
    var inputs_ : inputsConfig;
    var it_ : int = 0;

    var cells_dom_ : domain(1) = {1..0};

    var dtCells_ : [cells_dom_] real(64);

    var rhorho_0_ : [cells_dom_] real(64);

    var Rd0_0_ : [cells_dom_] real(64);

    var R0_ : [cells_dom_] real(64);

    const a1 = 0.25; const b1 = 1.0;
    const a2 = 0.1667; const b2 = 0.0;
    const a3 = 0.3750; const b3 = 0.56;
    const a4 = 0.5; const b4 = 0.0;
    const a5 = 1.0; const b5 = 0.44;

    var timeList = new list(real(64));
    var itList = new list(int);
    var res0List = new list(real(64));
    var clList = new list(real(64));
    var cdList = new list(real(64));
    var cmList = new list(real(64));

    proc init(spatialDisc : shared fullPotentialEikonalSpatialDiscretization, ref inputs : inputsConfig) {
        this.spatialDisc_ = spatialDisc;
        this.inputs_ = inputs;

        this.cells_dom_ = this.spatialDisc_.cell_dom_with_ghosts_;
    }

    proc compute_dt() {
        this.dtCells_ = 1.0e20;
        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            for i in 0..<this.spatialDisc_.mesh_.niCell_ {
                const cellIndexFVM = this.spatialDisc_.meshIndex2FVMindex(i, j);
                if (this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != 1) {
                    continue;
                }
                const lambdaI = this.spatialDisc_.LambdaI_[cellIndexFVM];
                const lambdaJ = this.spatialDisc_.LambdaJ_[cellIndexFVM];
                const volume = this.spatialDisc_.cellVolumesWithGhosts_[cellIndexFVM];
                this.dtCells_[cellIndexFVM] = this.inputs_.CFL_ * volume / (lambdaI + lambdaJ);

            }
        }
        // const min_dt = min reduce this.dtCells_;
        // writeln("Min dt: ", min_dt);
        // this.dtCells_ = min_dt;
    }

    proc eulerStep() {
        this.rhorho_0_ = this.spatialDisc_.rhorho_;

        this.spatialDisc_.run_odd();
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.rhorho_[cell] = this.rhorho_0_[cell] - this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc0_[cell] - this.spatialDisc_.Rd0_[cell]);
            this.R0_[cell] = this.spatialDisc_.rhorho_[cell] - this.rhorho_0_[cell];
        }
    }

    proc RKstep() {
        this.rhorho_0_ = this.spatialDisc_.rhorho_;

        // Stage 1
        this.spatialDisc_.run_odd();
        this.Rd0_0_ = this.spatialDisc_.Rd0_;

        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.rhorho_[cell] = this.rhorho_0_[cell] - a1 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc0_[cell] - this.Rd0_0_[cell]);
        }

        // Stage 2
        this.spatialDisc_.run_even();
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.rhorho_[cell] = this.rhorho_0_[cell] - a2 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc0_[cell] - this.Rd0_0_[cell]);
        }

        // Stage 3 - update dissipation
        this.spatialDisc_.run_odd();
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.Rd0_0_[cell] = b3 * this.spatialDisc_.Rd0_[cell] + (1 - b3) * this.Rd0_0_[cell];
        }
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.rhorho_[cell] = this.rhorho_0_[cell] - a3 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc0_[cell] - this.Rd0_0_[cell]);
        }

        // Stage 4
        this.spatialDisc_.run_even();
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.rhorho_[cell] = this.rhorho_0_[cell] - a4 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc0_[cell] - this.Rd0_0_[cell]);
        }

        // Stage 5 - update dissipation
        this.spatialDisc_.run_odd();
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.Rd0_0_[cell] = b5 * this.spatialDisc_.Rd0_[cell] + (1 - b5) * this.Rd0_0_[cell];
        }
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.rhorho_[cell] = this.rhorho_0_[cell] - a5 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc0_[cell] - this.Rd0_0_[cell]);

            this.R0_[cell] = this.spatialDisc_.rhorho_[cell] - this.rhorho_0_[cell];
        }
    }

    proc solve() {
        var time: stopwatch;
        var res0 = 1e12;
        var firstRes0 = 0.0;
        var firstRes1 = 0.0;
        var firstRes2 = 0.0;
        var firstRes3 = 0.0;

        while (res0 > this.inputs_.CONV_TOL_ && this.it_ < this.inputs_.IT_MAX_ && isNan(res0) == false) {
            this.it_ += 1;
            
            time.start();
            this.compute_dt();
            this.eulerStep();
            this.spatialDisc_.eikonalSolve();

            const (Cl, Cd, Cm) = this.spatialDisc_.compute_aerodynamics_coefficients();

            time.stop();

            res0 = l2Norm(this.R0_);

            if this.it_ == 1 {
                firstRes0 = res0;
            }

            res0 = res0 / firstRes0;

            writeln("t: ", time.elapsed()," Iteration: ", this.it_, " Residuals: ", res0, " Cl: ", Cl, " Cd: ", Cd, " Cm: ", Cm);

            this.timeList.pushBack(time.elapsed());
            this.itList.pushBack(this.it_);
            this.res0List.pushBack(res0);
            this.clList.pushBack(Cl);
            this.cdList.pushBack(Cd);
            this.cmList.pushBack(Cm);

            if this.it_ % this.inputs_.CGNS_OUTPUT_FREQ_ == 0 {
                writeln("Writing CGNS output at iteration ", this.it_);
                this.spatialDisc_.R0_ = this.R0_;
                this.spatialDisc_.writeSolution2CGNS(timeList, itList, res0List, clList, cdList, cmList);
            }
        }
        this.spatialDisc_.R0_ = this.R0_;
        this.spatialDisc_.writeSolution2CGNS(timeList, itList, res0List, clList, cdList, cmList);
    }
}

} // module temporalDiscretization