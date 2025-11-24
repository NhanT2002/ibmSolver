module fullPotentialTemporalDiscretization 
{
use mesh;
use fullPotentialSpatialDiscretization;
use linearAlgebra;
import input.inputsConfig;
use Time;
use Math;
use List;
use PETSCapi;
use C_PETSC;
use petsc;
use CTypes;

class fullPotentialTemporalDiscretization {
    var spatialDisc_ : shared fullPotentialSpatialDiscretization;
    var inputs_ : inputsConfig;
    var it_ : int = 0;

    var cells_dom_ : domain(1) = {1..0};

    var dtCells_ : [cells_dom_] real(64);

    var rho_0_ : [cells_dom_] real(64);
    var phi_0_ : [cells_dom_] real(64);

    var Rd0_0_ : [cells_dom_] real(64);
    var Rd1_0_ : [cells_dom_] real(64);

    var R0_ : [cells_dom_] real(64);
    var R1_ : [cells_dom_] real(64);

    const a1 = 0.25; const b1 = 1.0;
    const a2 = 0.1667; const b2 = 0.0;
    const a3 = 0.3750; const b3 = 0.56;
    const a4 = 0.5; const b4 = 0.0;
    const a5 = 1.0; const b5 = 0.44;

    var A_petsc : owned PETSCmatrix_c;
    var x_petsc : owned PETSCvector_c;
    var b_petsc : owned PETSCvector_c;
    var ksp : owned PETSCksp_c;

    var timeList = new list(real(64));
    var itList = new list(int);
    var res0List = new list(real(64));
    var res1List = new list(real(64));
    var res2List = new list(real(64));
    var res3List = new list(real(64));
    var clList = new list(real(64));
    var cdList = new list(real(64));
    var cmList = new list(real(64));

    proc init(spatialDisc : shared fullPotentialSpatialDiscretization, ref inputs : inputsConfig) {
        this.spatialDisc_ = spatialDisc;
        this.inputs_ = inputs;

        this.cells_dom_ = this.spatialDisc_.cell_dom_with_ghosts_;

        const M = this.spatialDisc_.mesh_.nCell_*2;
        const N = this.spatialDisc_.mesh_.nCell_*2;
        this.A_petsc = new owned PETSCmatrix_c(PETSC_COMM_SELF, "seqaij", M, M, N, N);
        this.x_petsc = new owned PETSCvector_c(PETSC_COMM_SELF, N, N, 0.0, "seq");
        this.b_petsc = new owned PETSCvector_c(PETSC_COMM_SELF, N, N, 0.0, "seq");

        var nnz : [0..M-1] PetscInt;
        nnz = 10; // Approximate number of non-zeros per row;
        A_petsc.preAllocate(nnz);

        this.ksp = new owned PETSCksp_c(PETSC_COMM_SELF, "gmres");
        this.ksp.setTolerances(inputs.GMRES_RTOL_, inputs.GMRES_ATOL_, inputs.GMRES_DTOL_, inputs.GMRES_MAXIT_);
        this.ksp.GMRESSetRestart(inputs.GMRES_RESTART_);
        this.ksp.GMRESSetPreAllocateVectors();
        if this.inputs_.GMRES_PRECON_ == "jacobi" {
            writeln("Using Jacobi preconditioner for GMRES");
            this.ksp.setPreconditioner("jacobi");
        } else if this.inputs_.GMRES_PRECON_ == "ilu" {
            writeln("Using ILU preconditioner for GMRES");
            this.ksp.setPreconditioner("ilu");
        } else if this.inputs_.GMRES_PRECON_ == "lu" {
            writeln("Using lu preconditioner for GMRES");
            this.ksp.setPreconditioner("lu");
        } else if this.inputs_.GMRES_PRECON_ == "asm" {
            writeln("Using asm preconditioner for GMRES");
            this.ksp.setPreconditioner("asm");
        } else if this.inputs_.GMRES_PRECON_ == "gasm" {
            writeln("Using gasm preconditioner for GMRES");
            this.ksp.setPreconditioner("gasm");
        } else if this.inputs_.GMRES_PRECON_ == "bjacobi" {
            writeln("Using bjacobi preconditioner for GMRES");
            this.ksp.setPreconditioner("bjacobi");
        } else if this.inputs_.GMRES_PRECON_ == "none" {
            writeln("Using no preconditioner for GMRES");
            this.ksp.setPreconditioner("none");
        } else {
            writeln("No preconditioner for GMRES");
            this.ksp.setPreconditioner("none");
        }
    }

    proc initializeJacobian() {
        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            for i in 0..<this.spatialDisc_.mesh_.niCell_ {
                const cellIndex = this.spatialDisc_.mesh_.iiCell(i, j);
                const cellIndexFVM = this.spatialDisc_.meshIndex2FVMindex(i, j);

                this.A_petsc.set(2*cellIndex, 2*cellIndex, 0.0);
                this.A_petsc.set(2*cellIndex, 2*cellIndex+1, 0.0);
                this.A_petsc.set(2*cellIndex+1, 2*cellIndex, 0.0);
                this.A_petsc.set(2*cellIndex+1, 2*cellIndex+1, 0.0);

                const leftCell = cellIndex - 1;
                const leftCellFVM = cellIndexFVM - 1;
                if this.spatialDisc_.cellTypesWithGhosts_[leftCellFVM] == 1  || this.spatialDisc_.cellTypesWithGhosts_[leftCellFVM] == 0 || this.spatialDisc_.cellTypesWithGhosts_[leftCellFVM] == -2 {
                    this.A_petsc.set(2*cellIndex, 2*leftCell, 0.0);
                    this.A_petsc.set(2*cellIndex, 2*leftCell+1, 0.0);
                    this.A_petsc.set(2*cellIndex+1, 2*leftCell+1, 0.0);
                }

                const rightCell = cellIndex + 1;
                const rightCellFVM = cellIndexFVM + 1;
                if this.spatialDisc_.cellTypesWithGhosts_[rightCellFVM] == 1 || this.spatialDisc_.cellTypesWithGhosts_[rightCellFVM] == 0 || this.spatialDisc_.cellTypesWithGhosts_[rightCellFVM] == -2 {
                    this.A_petsc.set(2*cellIndex, 2*rightCell, 0.0);
                    this.A_petsc.set(2*cellIndex, 2*rightCell+1, 0.0);
                    this.A_petsc.set(2*cellIndex+1, 2*rightCell+1, 0.0);
                }

                const bottomCell = cellIndex - this.spatialDisc_.mesh_.niCell_;
                const bottomCellFVM = cellIndexFVM - this.spatialDisc_.niCellWithGhosts_;
                if this.spatialDisc_.cellTypesWithGhosts_[bottomCellFVM] == 1 || this.spatialDisc_.cellTypesWithGhosts_[bottomCellFVM] == 0 || this.spatialDisc_.cellTypesWithGhosts_[bottomCellFVM] == -2 {
                    this.A_petsc.set(2*cellIndex, 2*bottomCell, 0.0);
                    this.A_petsc.set(2*cellIndex, 2*bottomCell+1, 0.0);
                    this.A_petsc.set(2*cellIndex+1, 2*bottomCell+1, 0.0);
                }

                const topCell = cellIndex + this.spatialDisc_.mesh_.niCell_;
                const topCellFVM = cellIndexFVM + this.spatialDisc_.niCellWithGhosts_;
                if this.spatialDisc_.cellTypesWithGhosts_[topCellFVM] == 1 || this.spatialDisc_.cellTypesWithGhosts_[topCellFVM] == 0 || this.spatialDisc_.cellTypesWithGhosts_[topCellFVM] == -2 {
                    this.A_petsc.set(2*cellIndex, 2*topCell, 0.0);
                    this.A_petsc.set(2*cellIndex, 2*topCell+1, 0.0);
                    this.A_petsc.set(2*cellIndex+1, 2*topCell+1, 0.0);
                }
            }
        }

        this.A_petsc.assemblyComplete();
        // this.A_petsc.matView();
    }

    proc computeJacobian() {
        this.A_petsc.zeroEntries();

        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            for i in 0..<this.spatialDisc_.mesh_.niCell_ {
                const cellIndex = this.spatialDisc_.mesh_.iiCell(i, j);
                const cellIndexFVM = this.spatialDisc_.meshIndex2FVMindex(i, j);
                const dt = this.dtCells_[cellIndexFVM];
                if (this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != 1 && this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != 0 && this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != -2) {
                    // Diagonal set to 1
                    this.A_petsc.set(2*cellIndex, 2*cellIndex, 1.0);
                    this.A_petsc.set(2*cellIndex+1, 2*cellIndex+1, 1.0);
                    continue;
                }

                const leftCellFVM = cellIndexFVM - 1;
                const rightCellFVM = cellIndexFVM + 1;
                const bottomCellFVM = cellIndexFVM - this.spatialDisc_.niCellWithGhosts_;
                const topCellFVM = cellIndexFVM + this.spatialDisc_.niCellWithGhosts_;

                const leftFace = i + j*this.spatialDisc_.mesh_.niNode_;
                const rightFace = leftFace + 1;
                const bottomFace = i + j*this.spatialDisc_.mesh_.niCell_;
                const topFace = bottomFace + this.spatialDisc_.mesh_.niCell_;

                // writeln("Cell (", i, ", ", j, ") Index: ", cellIndex, " FVM Index: ", cellIndexFVM, 
                // " Left : ", leftCellFVM, " Right : ", rightCellFVM, 
                // " Bottom : ", bottomCellFVM, " Top : ", topCellFVM);
                // writeln("leftFace: ", leftFace, " rightFace: ", rightFace, " bottomFace: ", bottomFace, " topFace: ", topFace);
                
                // Diagonal terms
                this.A_petsc.set(2*cellIndex, 2*cellIndex, this.spatialDisc_.cellVolumesWithGhosts_[cellIndexFVM]/dt);
                this.A_petsc.set(2*cellIndex+1, 2*cellIndex+1, this.spatialDisc_.cellVolumesWithGhosts_[cellIndexFVM]/dt);

                // dR0/drho_i terms
                this.A_petsc.add(2*cellIndex, 2*cellIndex, -this.spatialDisc_.dF0Idrho_[leftFace] + this.spatialDisc_.dF0Idrho_[rightFace] - this.spatialDisc_.dF0Jdrho_[bottomFace] + this.spatialDisc_.dF0Jdrho_[topFace]);

                // // dR0/drho_j terms
                // if this.spatialDisc_.cellTypesWithGhosts_[leftCellFVM] == 1 {
                //     const leftCell = cellIndex - 1;
                //     this.A_petsc.set(2*cellIndex, 2*leftCell, -this.spatialDisc_.dF0Idrho_[leftFace]);
                // }
                // if this.spatialDisc_.cellTypesWithGhosts_[rightCellFVM] == 1 {
                //     const rightCell = cellIndex + 1;
                //     this.A_petsc.set(2*cellIndex, 2*rightCell, this.spatialDisc_.dF0Idrho_[rightFace]);
                // }
                // if this.spatialDisc_.cellTypesWithGhosts_[bottomCellFVM] == 1 {
                //     const bottomCell = cellIndex - this.spatialDisc_.mesh_.niCell_;
                //     this.A_petsc.set(2*cellIndex, 2*bottomCell, -this.spatialDisc_.dF0Jdrho_[bottomFace]);
                // }
                // if this.spatialDisc_.cellTypesWithGhosts_[topCellFVM] == 1 {
                //     const topCell = cellIndex + this.spatialDisc_.mesh_.niCell_;
                //     this.A_petsc.set(2*cellIndex, 2*topCell, this.spatialDisc_.dF0Jdrho_[topFace]);
                // }

                // dR0/dphi_i terms
                const avg_rho_left = this.spatialDisc_.avg_rho_I_[leftFace];
                const avg_rho_right = this.spatialDisc_.avg_rho_I_[rightFace];
                const avg_rho_bottom = this.spatialDisc_.avg_rho_J_[bottomFace];
                const avg_rho_top = this.spatialDisc_.avg_rho_J_[topFace];

                const dF0dphi_left = avg_rho_left / (this.spatialDisc_.xCellsWithGhosts_[cellIndexFVM] - this.spatialDisc_.xCellsWithGhosts_[leftCellFVM]) * this.spatialDisc_.mesh_.IfaceAreas_[leftFace];
                const dF0dphi_right = - avg_rho_right / (this.spatialDisc_.xCellsWithGhosts_[rightCellFVM] - this.spatialDisc_.xCellsWithGhosts_[cellIndexFVM]) * this.spatialDisc_.mesh_.IfaceAreas_[rightFace];
                const dF0dphi_bottom = avg_rho_bottom / (this.spatialDisc_.yCellsWithGhosts_[cellIndexFVM] - this.spatialDisc_.yCellsWithGhosts_[bottomCellFVM]) * this.spatialDisc_.mesh_.JfaceAreas_[bottomFace];
                const dF0dphi_top = - avg_rho_top / (this.spatialDisc_.yCellsWithGhosts_[topCellFVM] - this.spatialDisc_.yCellsWithGhosts_[cellIndexFVM]) * this.spatialDisc_.mesh_.JfaceAreas_[topFace];

                this.A_petsc.add(2*cellIndex, 2*cellIndex+1, -dF0dphi_left + dF0dphi_right - dF0dphi_bottom + dF0dphi_top);

                // dR1/drho_i terms
                this.A_petsc.set(2*cellIndex+1, 2*cellIndex, this.spatialDisc_.cellVolumesWithGhosts_[cellIndexFVM]*this.spatialDisc_.rhorho_[cellIndexFVM]**(this.inputs_.GAMMA_-2.0) / this.inputs_.MACH_**2);

                // dR1/dphi_i terms
                this.A_petsc.add(2*cellIndex+1, 2*cellIndex+1, this.spatialDisc_.cellVolumesWithGhosts_[cellIndexFVM]*(this.spatialDisc_.uu_[cellIndexFVM] * this.spatialDisc_.si_[cellIndexFVM] + this.spatialDisc_.vv_[cellIndexFVM] * this.spatialDisc_.sj_[cellIndexFVM]));

                
                if this.spatialDisc_.cellTypesWithGhosts_[leftCellFVM] == 1 || this.spatialDisc_.cellTypesWithGhosts_[leftCellFVM] == 0 || this.spatialDisc_.cellTypesWithGhosts_[leftCellFVM] == -2 {
                    const leftCell = cellIndex - 1;
                    // dR1/dphi_j terms
                    this.A_petsc.set(2*cellIndex+1, 2*leftCell+1, this.spatialDisc_.cellVolumesWithGhosts_[cellIndexFVM]*(this.spatialDisc_.uu_[cellIndexFVM] * this.spatialDisc_.mi_[cellIndexFVM]));
                    // dR0/dphi_j terms
                    this.A_petsc.set(2*cellIndex, 2*leftCell+1, dF0dphi_left);
                }
                if this.spatialDisc_.cellTypesWithGhosts_[rightCellFVM] == 1 || this.spatialDisc_.cellTypesWithGhosts_[rightCellFVM] == 0 || this.spatialDisc_.cellTypesWithGhosts_[rightCellFVM] == -2 {
                    const rightCell = cellIndex + 1;
                    // dR1/dphi_j terms
                    this.A_petsc.set(2*cellIndex+1, 2*rightCell+1, this.spatialDisc_.cellVolumesWithGhosts_[cellIndexFVM]*(this.spatialDisc_.uu_[cellIndexFVM] * this.spatialDisc_.ni_[cellIndexFVM]));
                    // dR0/dphi_j terms
                    this.A_petsc.set(2*cellIndex, 2*rightCell+1, -dF0dphi_right);
                }
                if this.spatialDisc_.cellTypesWithGhosts_[bottomCellFVM] == 1 || this.spatialDisc_.cellTypesWithGhosts_[bottomCellFVM] == 0 || this.spatialDisc_.cellTypesWithGhosts_[bottomCellFVM] == -2 {
                    const bottomCell = cellIndex - this.spatialDisc_.mesh_.niCell_;
                    // dR1/dphi_j terms
                    this.A_petsc.set(2*cellIndex+1, 2*bottomCell+1, this.spatialDisc_.cellVolumesWithGhosts_[cellIndexFVM]*(this.spatialDisc_.vv_[cellIndexFVM] * this.spatialDisc_.mj_[cellIndexFVM]));
                    // dR0/dphi_j terms
                    this.A_petsc.set(2*cellIndex, 2*bottomCell+1, dF0dphi_bottom);
                }
                if this.spatialDisc_.cellTypesWithGhosts_[topCellFVM] == 1 || this.spatialDisc_.cellTypesWithGhosts_[topCellFVM] == 0 || this.spatialDisc_.cellTypesWithGhosts_[topCellFVM] == -2 {
                    const topCell = cellIndex + this.spatialDisc_.mesh_.niCell_;
                    // dR1/dphi_j terms
                    this.A_petsc.set(2*cellIndex+1, 2*topCell+1, this.spatialDisc_.cellVolumesWithGhosts_[cellIndexFVM]*(this.spatialDisc_.vv_[cellIndexFVM] * this.spatialDisc_.nj_[cellIndexFVM]));
                    // dR0/dphi_j terms
                    this.A_petsc.set(2*cellIndex, 2*topCell+1, -dF0dphi_top);
                }
            }
        }
        this.A_petsc.assemblyComplete();
        // this.A_petsc.matView();
    }

    proc compute_dt() {
        this.dtCells_ = this.inputs_.CFL_;
        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            for i in 0..<this.spatialDisc_.mesh_.niCell_ {
                const cellIndexFVM = this.spatialDisc_.meshIndex2FVMindex(i, j);
                if (this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != 1 && this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != 0 && this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != -2) {
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

        this.spatialDisc_.run();
        this.compute_dt();
        this.computeJacobian();
        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            for i in 0..<this.spatialDisc_.mesh_.niCell_ {
                const cellIndex = this.spatialDisc_.mesh_.iiCell(i, j);
                const cellIndexFVM = this.spatialDisc_.meshIndex2FVMindex(i, j);

                this.b_petsc.set(2*cellIndex, -this.spatialDisc_.Rc0_[cellIndexFVM]);
                this.b_petsc.set(2*cellIndex+1, -this.spatialDisc_.Rc1_[cellIndexFVM]);
            }
        }
        this.b_petsc.assemblyComplete();

        const (its, reason) = GMRES(this.ksp, this.A_petsc, this.b_petsc, this.x_petsc);

        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            for i in 0..<this.spatialDisc_.mesh_.niCell_ {
                const cellIndex = this.spatialDisc_.mesh_.iiCell(i, j);
                const cellIndexFVM = this.spatialDisc_.meshIndex2FVMindex(i, j);
                this.spatialDisc_.rhorho_[cellIndexFVM] += this.inputs_.OMEGA_ * this.x_petsc.get(2*cellIndex);
                this.spatialDisc_.phiphi_[cellIndexFVM] += this.inputs_.OMEGA_ * this.x_petsc.get(2*cellIndex+1);
                this.spatialDisc_.phi_t_[cellIndexFVM] = this.inputs_.OMEGA_ * this.x_petsc.get(2*cellIndex+1) / this.dtCells_[cellIndexFVM];

                // if cellIndex == 22015 {
                //     writeln("  cellIndex = ", cellIndex,
                //             " uu_ = ", this.spatialDisc_.uu_[cellIndexFVM],
                //             " vv_ = ", this.spatialDisc_.vv_[cellIndexFVM],
                //             " rho = ", this.spatialDisc_.rhorho_[cellIndexFVM],
                //             " phi = ", this.spatialDisc_.phiphi_[cellIndexFVM],
                //             " Rc0 = ", this.spatialDisc_.Rc0_[cellIndexFVM],
                //             " Rc1 = ", this.spatialDisc_.Rc1_[cellIndexFVM],
                //             " delta rho = ", this.inputs_.OMEGA_ * this.x_petsc.get(2*cellIndex),
                //             " delta phi = ", this.inputs_.OMEGA_ * this.x_petsc.get(2*cellIndex+1));
                //             writeln(" ");
                // }
            }
        }

        return (its, reason);
    }

    proc RKstep() {
        this.rho_0_ = this.spatialDisc_.rhorho_;
        this.phi_0_ = this.spatialDisc_.phiphi_;

        // Stage 1
        this.spatialDisc_.run();

        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.rhorho_[cell] = this.rho_0_[cell] - a1 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * this.spatialDisc_.Rc0_[cell];
            this.spatialDisc_.phiphi_[cell] = this.phi_0_[cell] - a1 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * this.spatialDisc_.Rc1_[cell];
        }

        // Stage 2
        this.spatialDisc_.run();
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.rhorho_[cell] = this.rho_0_[cell] - a2 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * this.spatialDisc_.Rc0_[cell];
            this.spatialDisc_.phiphi_[cell] = this.phi_0_[cell] - a2 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * this.spatialDisc_.Rc1_[cell];
        }

        // Stage 3 - update dissipation
        this.spatialDisc_.run();
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.rhorho_[cell] = this.rho_0_[cell] - a3 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * this.spatialDisc_.Rc0_[cell];
            this.spatialDisc_.phiphi_[cell] = this.phi_0_[cell] - a3 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * this.spatialDisc_.Rc1_[cell];
        }

        // Stage 4
        this.spatialDisc_.run();
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.rhorho_[cell] = this.rho_0_[cell] - a4 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * this.spatialDisc_.Rc0_[cell];
            this.spatialDisc_.phiphi_[cell] = this.phi_0_[cell] - a4 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * this.spatialDisc_.Rc1_[cell];
        }

        // Stage 5 - update dissipation
        this.spatialDisc_.run();
        forall cell in 0..<this.spatialDisc_.nCellWithGhosts_ {
            if (this.spatialDisc_.cellTypesWithGhosts_[cell] != 1) {
                continue;
            }
            this.spatialDisc_.rhorho_[cell] = this.rho_0_[cell] - a5 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * this.spatialDisc_.Rc0_[cell];
            this.spatialDisc_.phiphi_[cell] = this.phi_0_[cell] - a5 * this.dtCells_[cell] / this.spatialDisc_.cellVolumesWithGhosts_[cell] * this.spatialDisc_.Rc1_[cell];
        }
    }

    proc solve() {
        var time: stopwatch;
        var normalize_res0 = 1e12;
        var normalize_res1 = 1e12;
        var current_normalize_res0 = 1e12;
        var firstRes0 = 0.0;
        var firstRes1 = 0.0;

        this.initializeJacobian();
        
        while (normalize_res0 > this.inputs_.CONV_TOL_ && this.it_ < this.inputs_.IT_MAX_ && isNan(normalize_res0) == false && isNan(normalize_res1) == false) {
            this.it_ += 1;
            
            time.start();
            const (its, reason) = this.eulerStep();

            const (Cl, Cd, Cm) = this.spatialDisc_.compute_aerodynamics_coefficients();

            time.stop();

            const res0 = l2Norm(this.spatialDisc_.Rc0_);
            const res1 = l2Norm(this.spatialDisc_.Rc1_);
            const dphidt_norm = l2Norm(this.spatialDisc_.phi_t_);

            if this.it_ == 1 {
                firstRes0 = res0;
                firstRes1 = res1;
            }

            current_normalize_res0 = res0 / firstRes0;
            
            if this.it_ > this.inputs_.CFL_RAMP_IT_ {
                this.inputs_.CFL_ = min(this.inputs_.CFL_ * abs(normalize_res0 / current_normalize_res0), this.inputs_.CFL_RAMP_MAX_);
            }


            normalize_res0 = res0 / firstRes0;
            normalize_res1 = res1 / firstRes1;

            // if this.it_ % this.inputs_.CFL_RAMP_IT_ == 0 {
            //     this.inputs_.CFL_ = min(this.inputs_.CFL_ * this.inputs_.CFL_RAMP_FACTOR_, this.inputs_.CFL_RAMP_FINAL_);
            //     writeln("Ramping CFL to ", this.inputs_.CFL_, " at iteration ", this.it_);
            // }

            writeln("t: ", time.elapsed()," Iteration: ", this.it_, " Res: ", res0, ", ", res1, 
            " Cl: ", Cl, " Cd: ", Cd, " Cm: ", Cm, " GMRES its: ", its, " Reason: ", reason, 
            " Norm Res: ", normalize_res0, ", ", normalize_res1, " CFL: ", this.inputs_.CFL_, " dphidt norm: ", dphidt_norm);

            this.timeList.pushBack(time.elapsed());
            this.itList.pushBack(this.it_);
            this.res0List.pushBack(res0);
            this.res1List.pushBack(res1);
            this.clList.pushBack(Cl);
            this.cdList.pushBack(Cd);
            this.cmList.pushBack(Cm);

            if this.it_ % this.inputs_.CGNS_OUTPUT_FREQ_ == 0 {
                writeln("Writing CGNS output at iteration ", this.it_);
                this.spatialDisc_.writeSolution2CGNS(timeList, itList, res0List, res1List, clList, cdList, cmList);
            }
        }
        writeln("jacobian = ");
        this.A_petsc.matView();
        // writeln("final delta rho and delta phi ");
        // for j in 0..<this.spatialDisc_.mesh_.njCell_ {
        //     for i in 0..<this.spatialDisc_.mesh_.niCell_ {
        //         const cellIndex = this.spatialDisc_.mesh_.iiCell(i, j);
        //         const cellIndexFVM = this.spatialDisc_.meshIndex2FVMindex(i, j);

        //         writeln("  cellIndex = ", cellIndex,
        //                 " uu_ = ", this.spatialDisc_.uu_[cellIndexFVM],
        //                 " vv_ = ", this.spatialDisc_.vv_[cellIndexFVM],
        //                 " rho = ", this.spatialDisc_.rhorho_[cellIndexFVM],
        //                 " phi = ", this.spatialDisc_.phiphi_[cellIndexFVM],
        //                 " Rc0 = ", this.spatialDisc_.Rc0_[cellIndexFVM],
        //                 " Rc1 = ", this.spatialDisc_.Rc1_[cellIndexFVM],
        //                 " delta rho = ", this.inputs_.OMEGA_ * this.x_petsc.get(2*cellIndex),
        //                 " delta phi = ", this.inputs_.OMEGA_ * this.x_petsc.get(2*cellIndex+1));

        //     }
        // }
        this.spatialDisc_.writeSolution2CGNS(timeList, itList, res0List, res1List, clList, cdList, cmList);
    }
}

} // module fullPotentialTemporalDiscretization