module temporalDiscretization 
{
use mesh;
use spatialDiscretization;
use linearAlgebra;
import input.inputsConfig;
use Time;
use Math;
use List;

class temporalDiscretization {
    var spatialDisc_ : shared spatialDiscretization;
    var inputs_ : inputsConfig;
    var it_ : int = 0;

    var cells_dom_ : domain(1) = {1..0};

    var dtCells_ : [cells_dom_] real(64);

    var W0_0_ : [cells_dom_] real(64);
    var W1_0_ : [cells_dom_] real(64);
    var W2_0_ : [cells_dom_] real(64);
    var W3_0_ : [cells_dom_] real(64);

    var Rd0_0_ : [cells_dom_] real(64);
    var Rd1_0_ : [cells_dom_] real(64);
    var Rd2_0_ : [cells_dom_] real(64);
    var Rd3_0_ : [cells_dom_] real(64);

    var R0_ : [cells_dom_] real(64);
    var R1_ : [cells_dom_] real(64);
    var R2_ : [cells_dom_] real(64);
    var R3_ : [cells_dom_] real(64);

    const a1 = 0.25; const b1 = 1.0;
    const a2 = 0.1667; const b2 = 0.0;
    const a3 = 0.3750; const b3 = 0.56;
    const a4 = 0.5; const b4 = 0.0;
    const a5 = 1.0; const b5 = 0.44;

    var timeList = new list(real(64));
    var itList = new list(int);
    var res0List = new list(real(64));
    var res1List = new list(real(64));
    var res2List = new list(real(64));
    var res3List = new list(real(64));
    var clList = new list(real(64));
    var cdList = new list(real(64));
    var cmList = new list(real(64));

    proc init(spatialDisc : shared spatialDiscretization, ref inputs : inputsConfig) {
        this.spatialDisc_ = spatialDisc;
        this.inputs_ = inputs;

        this.cells_dom_ = this.spatialDisc_.cell_dom_with_ghosts_;
    }

    proc compute_dt() {
        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            for i in 0..<this.spatialDisc_.mesh_.niCell_ {
                const cellIndexFVM = this.spatialDisc_.meshIndex2FVMindex(i, j);
                const lambdaI = this.spatialDisc_.LambdaI_[cellIndexFVM];
                const lambdaJ = this.spatialDisc_.LambdaJ_[cellIndexFVM];
                const volume = this.spatialDisc_.cellVolumesWithGhosts_[cellIndexFVM];
                this.dtCells_[cellIndexFVM] = this.inputs_.CFL_ * volume / (lambdaI + lambdaJ);

            }
        }
    }

    proc RKstep() {
        this.W0_0_ = this.spatialDisc_.W0_;
        this.W1_0_ = this.spatialDisc_.W1_;
        this.W2_0_ = this.spatialDisc_.W2_;
        this.W3_0_ = this.spatialDisc_.W3_;

        // Stage 1
        this.spatialDisc_.run_odd();
        this.Rd0_0_ = this.spatialDisc_.Rd0_;
        this.Rd1_0_ = this.spatialDisc_.Rd1_;
        this.Rd2_0_ = this.spatialDisc_.Rd2_;
        this.Rd3_0_ = this.spatialDisc_.Rd3_;

        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.W0_[cell] = this.W0_0_[cell] - a1 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc0_[cell] - this.Rd0_0_[cell]);
            this.spatialDisc_.W1_[cell] = this.W1_0_[cell] - a1 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc1_[cell] - this.Rd1_0_[cell]);
            this.spatialDisc_.W2_[cell] = this.W2_0_[cell] - a1 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc2_[cell] - this.Rd2_0_[cell]);
            this.spatialDisc_.W3_[cell] = this.W3_0_[cell] - a1 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc3_[cell] - this.Rd3_0_[cell]);
        }

        // Stage 2
        this.spatialDisc_.run_even();
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.W0_[cell] = this.W0_0_[cell] - a2 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc0_[cell] - this.Rd0_0_[cell]);
            this.spatialDisc_.W1_[cell] = this.W1_0_[cell] - a2 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc1_[cell] - this.Rd1_0_[cell]);
            this.spatialDisc_.W2_[cell] = this.W2_0_[cell] - a2 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc2_[cell] - this.Rd2_0_[cell]);
            this.spatialDisc_.W3_[cell] = this.W3_0_[cell] - a2 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc3_[cell] - this.Rd3_0_[cell]);
        }

        // Stage 3 - update dissipation
        this.spatialDisc_.run_odd();
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.Rd0_0_[cell] = b3 * this.spatialDisc_.Rd0_[cell] + (1 - b3) * this.Rd0_0_[cell];
            this.Rd1_0_[cell] = b3 * this.spatialDisc_.Rd1_[cell] + (1 - b3) * this.Rd1_0_[cell];
            this.Rd2_0_[cell] = b3 * this.spatialDisc_.Rd2_[cell] + (1 - b3) * this.Rd2_0_[cell];
            this.Rd3_0_[cell] = b3 * this.spatialDisc_.Rd3_[cell] + (1 - b3) * this.Rd3_0_[cell];
        }
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.W0_[cell] = this.W0_0_[cell] - a3 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc0_[cell] - this.Rd0_0_[cell]);
            this.spatialDisc_.W1_[cell] = this.W1_0_[cell] - a3 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc1_[cell] - this.Rd1_0_[cell]);
            this.spatialDisc_.W2_[cell] = this.W2_0_[cell] - a3 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc2_[cell] - this.Rd2_0_[cell]);
            this.spatialDisc_.W3_[cell] = this.W3_0_[cell] - a3 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc3_[cell] - this.Rd3_0_[cell]);
        }

        // Stage 4
        this.spatialDisc_.run_even();
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.W0_[cell] = this.W0_0_[cell] - a4 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc0_[cell] - this.Rd0_0_[cell]);
            this.spatialDisc_.W1_[cell] = this.W1_0_[cell] - a4 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc1_[cell] - this.Rd1_0_[cell]);
            this.spatialDisc_.W2_[cell] = this.W2_0_[cell] - a4 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc2_[cell] - this.Rd2_0_[cell]);
            this.spatialDisc_.W3_[cell] = this.W3_0_[cell] - a4 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc3_[cell] - this.Rd3_0_[cell]);
        }

        // Stage 5 - update dissipation
        this.spatialDisc_.run_odd();
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.Rd0_0_[cell] = b5 * this.spatialDisc_.Rd0_[cell] + (1 - b5) * this.Rd0_0_[cell];
            this.Rd1_0_[cell] = b5 * this.spatialDisc_.Rd1_[cell] + (1 - b5) * this.Rd1_0_[cell];
            this.Rd2_0_[cell] = b5 * this.spatialDisc_.Rd2_[cell] + (1 - b5) * this.Rd2_0_[cell];
            this.Rd3_0_[cell] = b5 * this.spatialDisc_.Rd3_[cell] + (1 - b5) * this.Rd3_0_[cell];
        }
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.W0_[cell] = this.W0_0_[cell] - a5 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc0_[cell] - this.Rd0_0_[cell]);
            this.spatialDisc_.W1_[cell] = this.W1_0_[cell] - a5 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc1_[cell] - this.Rd1_0_[cell]);
            this.spatialDisc_.W2_[cell] = this.W2_0_[cell] - a5 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc2_[cell] - this.Rd2_0_[cell]);
            this.spatialDisc_.W3_[cell] = this.W3_0_[cell] - a5 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * (this.spatialDisc_.Rc3_[cell] - this.Rd3_0_[cell]);

            this.R0_[cell] = this.spatialDisc_.Rc0_[cell] - this.Rd0_0_[cell];
            this.R1_[cell] = this.spatialDisc_.Rc1_[cell] - this.Rd1_0_[cell];
            this.R2_[cell] = this.spatialDisc_.Rc2_[cell] - this.Rd2_0_[cell];
            this.R3_[cell] = this.spatialDisc_.Rc3_[cell] - this.Rd3_0_[cell];
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
            this.RKstep();

            const (Cl, Cd, Cm) = this.spatialDisc_.compute_aerodynamics_coefficients();

            time.stop();

            res0 = l2Norm(this.R0_);
            var res1 = l2Norm(this.R1_);
            var res2 = l2Norm(this.R2_);
            var res3 = l2Norm(this.R3_);

            if this.it_ == 1 {
                firstRes0 = res0;
                firstRes1 = res1;
                firstRes2 = res2;
                firstRes3 = res3;
            }

            res0 = res0 / firstRes0;
            res1 = res1 / firstRes1;
            res2 = res2 / firstRes2;
            res3 = res3 / firstRes3;

            writeln("t: ", time.elapsed()," Iteration: ", this.it_, " Residuals: ", res0, ", ", res1, ", ", res2, ", ", res3, " Cl: ", Cl, " Cd: ", Cd, " Cm: ", Cm);

            this.timeList.pushBack(time.elapsed());
            this.itList.pushBack(this.it_);
            this.res0List.pushBack(res0);
            this.res1List.pushBack(res1);
            this.res2List.pushBack(res2);
            this.res3List.pushBack(res3);
            this.clList.pushBack(Cl);
            this.cdList.pushBack(Cd);
            this.cmList.pushBack(Cm);

            if this.it_ % this.inputs_.CGNS_OUTPUT_FREQ_ == 0 {
                this.spatialDisc_.writeSolution2CGNS(timeList, itList, res0List, res1List, res2List, res3List, clList, cdList, cmList);
            }
        }
        this.spatialDisc_.writeSolution2CGNS(timeList, itList, res0List, res1List, res2List, res3List, clList, cdList, cmList);
    }
}

} // module temporalDiscretization