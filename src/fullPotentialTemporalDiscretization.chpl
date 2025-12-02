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
    var firstRes0 : real(64) = 0.0;

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

    var dom_ : domain(1) = {1..0};

    var Lx_a_ : [dom_] real(64);
    var Lx_b_ : [dom_] real(64);
    var Lx_c_ : [dom_] real(64);
    var Lx_d_ : [dom_] real(64);

    var Ly_a_ : [dom_] real(64);
    var Ly_b_ : [dom_] real(64);
    var Ly_c_ : [dom_] real(64);
    var Ly_d_ : [dom_] real(64);

    var x_ : [dom_] real(64);

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

    var wakeFaces_dom_ : domain(1) = {1..0};
    var deltaGamma_ : [wakeFaces_dom_] real(64);
    var wakeVelocities_ : [wakeFaces_dom_] real(64);
    var wakeCirculation_ : [wakeFaces_dom_] real(64);
    var oldCirculation_ = 0.0;

    var t_ : real(64) = 0.0;

    proc init(spatialDisc : shared fullPotentialSpatialDiscretization, ref inputs : inputsConfig) {
        this.spatialDisc_ = spatialDisc;
        this.inputs_ = inputs;

        this.cells_dom_ = this.spatialDisc_.cell_dom_with_ghosts_;
        this.dom_ = this.spatialDisc_.mesh_.cell_dom_;

        const M = this.spatialDisc_.mesh_.nCell_;
        const N = this.spatialDisc_.mesh_.nCell_;
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
                if (this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != 1 && this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != 0) {
                    // Diagonal set to 1
                    this.A_petsc.set(cellIndex, cellIndex, 1.0);
                    continue;
                }

                const leftCell = cellIndex - 1;
                const rightCell = cellIndex + 1;
                const bottomCell = cellIndex - this.spatialDisc_.mesh_.niCell_;
                const topCell = cellIndex + this.spatialDisc_.mesh_.niCell_;

                const leftCellFVM = cellIndexFVM - 1;
                const rightCellFVM = cellIndexFVM + 1;
                const bottomCellFVM = cellIndexFVM - this.spatialDisc_.niCellWithGhosts_;
                const topCellFVM = cellIndexFVM + this.spatialDisc_.niCellWithGhosts_;

                // diagonal term
                this.A_petsc.set(cellIndex, cellIndex, 0.0);

                // left neighbor
                if this.spatialDisc_.cellTypesWithGhosts_[leftCellFVM] != 9 {
                    this.A_petsc.set(cellIndex, leftCell, 0.0);
                }

                // right neighbor
                if this.spatialDisc_.cellTypesWithGhosts_[rightCellFVM] != 9 {
                    this.A_petsc.set(cellIndex, rightCell, 0.0);
                }

                // bottom neighbor
                if this.spatialDisc_.cellTypesWithGhosts_[bottomCellFVM] != 9 {
                    this.A_petsc.set(cellIndex, bottomCell, 0.0);
                }

                // top neighbor
                if this.spatialDisc_.cellTypesWithGhosts_[topCellFVM] != 9 {
                    this.A_petsc.set(cellIndex, topCell, 0.0);
                }
            }
        }

        this.A_petsc.assemblyComplete();
        // this.A_petsc.matView();
    }

    proc computeJacobian() {
        this.A_petsc.zeroEntries();
        const dt = this.inputs_.CFL_;

        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            for i in 0..<this.spatialDisc_.mesh_.niCell_ {
                const cellIndex = this.spatialDisc_.mesh_.iiCell(i, j);
                const cellIndexFVM = this.spatialDisc_.meshIndex2FVMindex(i, j);
                if (this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != 1 && this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != 0) {
                    // Diagonal set to 1
                    this.A_petsc.set(cellIndex, cellIndex, 1.0);
                    continue;
                }

                const leftCell = cellIndex - 1;
                const rightCell = cellIndex + 1;
                const bottomCell = cellIndex - this.spatialDisc_.mesh_.niCell_;
                const topCell = cellIndex + this.spatialDisc_.mesh_.niCell_;

                const leftCellFVM = cellIndexFVM - 1;
                const rightCellFVM = cellIndexFVM + 1;
                const bottomCellFVM = cellIndexFVM - this.spatialDisc_.niCellWithGhosts_;
                const topCellFVM = cellIndexFVM + this.spatialDisc_.niCellWithGhosts_;

                const u = this.spatialDisc_.uu_[cellIndexFVM];
                const v = this.spatialDisc_.vv_[cellIndexFVM];
                const rho = this.spatialDisc_.rhorho_[cellIndexFVM];
                const m_i = this.spatialDisc_.mi_[cellIndexFVM];
                const n_i = this.spatialDisc_.ni_[cellIndexFVM];
                const s_i = this.spatialDisc_.si_[cellIndexFVM];
                const m_j = this.spatialDisc_.mj_[cellIndexFVM];
                const n_j = this.spatialDisc_.nj_[cellIndexFVM];
                const s_j = this.spatialDisc_.sj_[cellIndexFVM];
                const beta = this.spatialDisc_.beta_[cellIndexFVM];

                const rho_left_face = 0.5*(this.spatialDisc_.rhorho_[leftCellFVM] + this.spatialDisc_.rhorho_[cellIndexFVM]);
                const rho_right_face = 0.5*(this.spatialDisc_.rhorho_[rightCellFVM] + this.spatialDisc_.rhorho_[cellIndexFVM]);
                const rho_bottom_face = 0.5*(this.spatialDisc_.rhorho_[bottomCellFVM] + this.spatialDisc_.rhorho_[cellIndexFVM]);
                const rho_top_face = 0.5*(this.spatialDisc_.rhorho_[topCellFVM] + this.spatialDisc_.rhorho_[cellIndexFVM]);

                const dx_left = this.spatialDisc_.xCellsWithGhosts_[cellIndexFVM] - this.spatialDisc_.xCellsWithGhosts_[leftCellFVM];
                const dx_right = this.spatialDisc_.xCellsWithGhosts_[rightCellFVM] - this.spatialDisc_.xCellsWithGhosts_[cellIndexFVM];
                const dy_bottom = this.spatialDisc_.yCellsWithGhosts_[cellIndexFVM] - this.spatialDisc_.yCellsWithGhosts_[bottomCellFVM];
                const dy_top = this.spatialDisc_.yCellsWithGhosts_[topCellFVM] - this.spatialDisc_.yCellsWithGhosts_[cellIndexFVM];

                const avg_dx_cell = this.spatialDisc_.mesh_.avgFaceAreaJ_[cellIndex];
                const avg_dy_cell = this.spatialDisc_.mesh_.avgFaceAreaI_[cellIndex];

                // diagonal term
                const Jij = 1.0 / dt**2 + (u*s_i + v*s_j) / dt - ((-rho_right_face/dx_right - rho_left_face/dx_left) / avg_dx_cell 
                                                                    + (-rho_top_face/dy_top - rho_bottom_face/dy_bottom) / avg_dy_cell) / beta;
                this.A_petsc.set(cellIndex, cellIndex, Jij);

                // left neighbor
                if this.spatialDisc_.cellTypesWithGhosts_[leftCellFVM] != 9 {
                    const Jim1j = (u*m_i) / dt - (rho_left_face/dx_left) / (beta * avg_dx_cell);
                    this.A_petsc.set(cellIndex, leftCell, Jim1j);
                }

                // right neighbor
                if this.spatialDisc_.cellTypesWithGhosts_[rightCellFVM] != 9 {
                    const Jip1j = (u*n_i) / dt - (rho_right_face/dx_right) / (beta * avg_dx_cell);
                    this.A_petsc.set(cellIndex, rightCell, Jip1j);
                }

                // bottom neighbor
                if this.spatialDisc_.cellTypesWithGhosts_[bottomCellFVM] != 9 {
                    const Jijm1 = (v*m_j) / dt - (rho_bottom_face/dy_bottom) / (beta * avg_dy_cell);
                    this.A_petsc.set(cellIndex, bottomCell, Jijm1);
                }

                // top neighbor
                if this.spatialDisc_.cellTypesWithGhosts_[topCellFVM] != 9 {
                    const Jijp1 = (v*n_j) / dt - (rho_top_face/dy_top) / (beta * avg_dy_cell);
                    this.A_petsc.set(cellIndex, topCell, Jijp1);
                }
            }
        }
        this.A_petsc.assemblyComplete();

        
        // forall cellIndex in 0..this.spatialDisc_.mesh_.nCell_-1 {
        //     if cellIndex == 25189 || cellIndex == 24965 {
        //         const leftCell = cellIndex - 1;
        //         const rightCell = cellIndex + 1;
        //         const bottomCell = cellIndex - this.spatialDisc_.mesh_.niCell_;;
        //         const topCell = cellIndex + this.spatialDisc_.mesh_.niCell_;
        //         writeln("cellIndex = ", cellIndex,
        //         " A[", cellIndex, ", ", cellIndex, "] = ", this.A_petsc.get(cellIndex, cellIndex),
        //         " A[", cellIndex, ", ", leftCell, "] = ", this.A_petsc.get(cellIndex, leftCell),
        //         " A[", cellIndex, ", ", rightCell, "] = ", this.A_petsc.get(cellIndex, rightCell),
        //         " A[", cellIndex, ", ", bottomCell, "] = ", this.A_petsc.get(cellIndex, bottomCell),
        //         " A[", cellIndex, ", ", topCell, "] = ", this.A_petsc.get(cellIndex, topCell));
        //     }
        // }
        // this.A_petsc.matView();
    }

    proc eulerStep() {

        this.spatialDisc_.run();
        this.computeJacobian();
        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            for i in 0..<this.spatialDisc_.mesh_.niCell_ {
                const cellIndex = this.spatialDisc_.mesh_.iiCell(i, j);
                const cellIndexFVM = this.spatialDisc_.meshIndex2FVMindex(i, j);

                this.b_petsc.set(cellIndex, this.spatialDisc_.R0_[cellIndexFVM]);
            }
        }
        this.b_petsc.assemblyComplete();

        const (its, reason) = GMRES(this.ksp, this.A_petsc, this.b_petsc, this.x_petsc);

        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            for i in 0..<this.spatialDisc_.mesh_.niCell_ {
                const cellIndex = this.spatialDisc_.mesh_.iiCell(i, j);
                const cellIndexFVM = this.spatialDisc_.meshIndex2FVMindex(i, j);

                if this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != 1 && this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != 0 {
                    continue;
                }
                
                // update past time level variables
                this.spatialDisc_.phiphi_m2_[cellIndexFVM] = this.spatialDisc_.phiphi_m1_[cellIndexFVM];
                
                this.spatialDisc_.rhorho_m1_[cellIndexFVM] = this.spatialDisc_.rhorho_[cellIndexFVM];
                this.spatialDisc_.beta_m1_[cellIndexFVM] = this.spatialDisc_.beta_[cellIndexFVM];
                this.spatialDisc_.phiphi_m1_[cellIndexFVM] = this.spatialDisc_.phiphi_[cellIndexFVM];
                this.spatialDisc_.uu_m1_[cellIndexFVM] = this.spatialDisc_.uu_[cellIndexFVM];
                this.spatialDisc_.vv_m1_[cellIndexFVM] = this.spatialDisc_.vv_[cellIndexFVM];

                // update solution
                this.spatialDisc_.phiphi_[cellIndexFVM] += this.inputs_.OMEGA_ * this.x_petsc.get(cellIndex);
                this.spatialDisc_.phi_minus_phi_m1_[cellIndexFVM] = this.spatialDisc_.phiphi_[cellIndexFVM] - this.spatialDisc_.phiphi_m1_[cellIndexFVM];

                // if cellIndex == 27756 {
                //     writeln("Updating cell ", cellIndex, ": phi = ", this.spatialDisc_.phiphi_[cellIndexFVM],
                //     " uu = ", this.spatialDisc_.uu_[cellIndexFVM],
                //     " vv = ", this.spatialDisc_.vv_[cellIndexFVM],
                //     " rho = ", this.spatialDisc_.rhorho_[cellIndexFVM],
                //     " delta phi = ", this.x_petsc.get(cellIndex),
                //     " residual = ", this.spatialDisc_.R0_[cellIndexFVM]);
                // }

            }
        }

        return (its, reason);
    }

    proc computeLx_Ly() {
        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            for i in 0..<this.spatialDisc_.mesh_.niCell_ {
                const cellIndex = this.spatialDisc_.mesh_.iiCell(i, j);
                const cellIndexFVM = this.spatialDisc_.meshIndex2FVMindex(i, j);
                if (this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != 1) {
                    // Diagonal set to 1
                    this.Lx_b_[cellIndex] = 1.0;
                    continue;
                }

                const dt = this.spatialDisc_.dtCells_[cellIndexFVM];
                const leftCellFVM = cellIndexFVM - 1;
                const rightCellFVM = cellIndexFVM + 1;

                const u = this.spatialDisc_.uu_[cellIndexFVM];
                const rho = this.spatialDisc_.rhorho_[cellIndexFVM];
                const m_i = this.spatialDisc_.mi_[cellIndexFVM];
                const n_i = this.spatialDisc_.ni_[cellIndexFVM];
                const s_i = this.spatialDisc_.si_[cellIndexFVM];
                const beta = this.spatialDisc_.beta_[cellIndexFVM];

                const rho_left_face = 0.5*(this.spatialDisc_.rhorho_[leftCellFVM] + this.spatialDisc_.rhorho_[cellIndexFVM]);
                const rho_right_face = 0.5*(this.spatialDisc_.rhorho_[rightCellFVM] + this.spatialDisc_.rhorho_[cellIndexFVM]);

                const dx_left = this.spatialDisc_.xCellsWithGhosts_[cellIndexFVM] - this.spatialDisc_.xCellsWithGhosts_[leftCellFVM];
                const dx_right = this.spatialDisc_.xCellsWithGhosts_[rightCellFVM] - this.spatialDisc_.xCellsWithGhosts_[cellIndexFVM];
                const avg_dx_cell = this.spatialDisc_.mesh_.avgFaceAreaJ_[cellIndex];

                // diagonal term
                const Jij_x = 1.0 / dt**2 + (u*s_i) / dt - ((-rho_right_face/dx_right - rho_left_face/dx_left) / avg_dx_cell) / beta;
                this.Lx_b_[cellIndex] = Jij_x;

                // left neighbor
                if this.spatialDisc_.cellTypesWithGhosts_[leftCellFVM] != 9 {
                    const Jim1j = (u*m_i) / dt - (rho_left_face/dx_left) / (beta * avg_dx_cell);
                    this.Lx_a_[cellIndex] = Jim1j;
                }

                // right neighbor
                if this.spatialDisc_.cellTypesWithGhosts_[rightCellFVM] != 9 {
                    const Jip1j = (u*n_i) / dt - (rho_right_face/dx_right) / (beta * avg_dx_cell);
                    this.Lx_c_[cellIndex] = Jip1j;
                }
            }
        }
        forall i in 0..<this.spatialDisc_.mesh_.niCell_ {
            for j in 0..<this.spatialDisc_.mesh_.njCell_ {
                const cellIndex = this.spatialDisc_.mesh_.iiCell(j, i);
                const cellIndexFVM = this.spatialDisc_.meshIndex2FVMindex(i, j);
                if (this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != 1) {
                    // Diagonal set to 1
                    this.Ly_b_[cellIndex] = 1.0;
                    continue;
                }

                const dt = this.spatialDisc_.dtCells_[cellIndexFVM];
                const bottomCellFVM = cellIndexFVM - this.spatialDisc_.niCellWithGhosts_;
                const topCellFVM = cellIndexFVM + this.spatialDisc_.niCellWithGhosts_;

                const v = this.spatialDisc_.vv_[cellIndexFVM];
                const rho = this.spatialDisc_.rhorho_[cellIndexFVM];
                const m_j = this.spatialDisc_.mj_[cellIndexFVM];
                const n_j = this.spatialDisc_.nj_[cellIndexFVM];
                const s_j = this.spatialDisc_.sj_[cellIndexFVM];
                const beta = this.spatialDisc_.beta_[cellIndexFVM];

                const rho_bottom_face = 0.5*(this.spatialDisc_.rhorho_[bottomCellFVM] + this.spatialDisc_.rhorho_[cellIndexFVM]);
                const rho_top_face = 0.5*(this.spatialDisc_.rhorho_[topCellFVM] + this.spatialDisc_.rhorho_[cellIndexFVM]);

                const dy_bottom = this.spatialDisc_.yCellsWithGhosts_[cellIndexFVM] - this.spatialDisc_.yCellsWithGhosts_[bottomCellFVM];
                const dy_top = this.spatialDisc_.yCellsWithGhosts_[topCellFVM] - this.spatialDisc_.yCellsWithGhosts_[cellIndexFVM];

                const avg_dy_cell = this.spatialDisc_.mesh_.avgFaceAreaI_[this.spatialDisc_.mesh_.iiCell(i, j)];

                // diagonal term
                const Jij_y = 1.0 / dt**2 + (v*s_j) / dt - ((-rho_top_face/dy_top - rho_bottom_face/dy_bottom) / avg_dy_cell) / beta;
                this.Ly_b_[cellIndex] = Jij_y;

                // bottom neighbor
                if this.spatialDisc_.cellTypesWithGhosts_[bottomCellFVM] != 9 {
                    const Jijm1 = (v*m_j) / dt - (rho_bottom_face/dy_bottom) / (beta * avg_dy_cell);
                    this.Ly_a_[cellIndex] = Jijm1;
                }

                // top neighbor
                if this.spatialDisc_.cellTypesWithGhosts_[topCellFVM] != 9 {
                    const Jijp1 = (v*n_j) / dt - (rho_top_face/dy_top) / (beta * avg_dy_cell);
                    this.Ly_c_[cellIndex] = Jijp1;
                }

            }
        }
    }

    proc printArray(inputArray : [dom_] real(64), niCell : int, njCell : int, arrayName : string) {
        writeln(arrayName, ":");
        for j in 0..<njCell {
            write("Row ", j, ": ");
            writeln(inputArray[this.spatialDisc_.mesh_.iiCell(0, j)..this.spatialDisc_.mesh_.iiCell(niCell, j)]);
        }
    }

    proc reArrangeRow2ColumnMajor(inputArray : [dom_] real(64)) {
        var outputArray : [dom_] real(64);

        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            for i in 0..<this.spatialDisc_.mesh_.niCell_ {
                const cellIndex = this.spatialDisc_.mesh_.iiCell(i, j);
                const cellIndexColumnMajor = this.spatialDisc_.mesh_.iiCell(j, i);

                outputArray[cellIndexColumnMajor] = inputArray[cellIndex];
            }
        }

        // printArray(inputArray, this.spatialDisc_.mesh_.niCell_, this.spatialDisc_.mesh_.njCell_, "inputArray (row-major)");
        // printArray(outputArray, this.spatialDisc_.mesh_.niCell_, this.spatialDisc_.mesh_.njCell_, "outputArray (column-major)");

        return outputArray;
    }

    proc eulerStepThomas() {
        this.spatialDisc_.run();
        this.computeLx_Ly();

        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            for i in 0..<this.spatialDisc_.mesh_.niCell_ {
                const cellIndex = this.spatialDisc_.mesh_.iiCell(i, j);
                const cellIndexFVM = this.spatialDisc_.meshIndex2FVMindex(i, j);

                this.Lx_d_[cellIndex] = this.spatialDisc_.R0_[cellIndexFVM];
            }
        }
        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            const cells = this.spatialDisc_.mesh_.iiCell(0, j)..this.spatialDisc_.mesh_.iiCell(this.spatialDisc_.mesh_.niCell_-1, j);
            thomasAlgorithm(this.Lx_a_[cells], this.Lx_b_[cells], this.Lx_c_[cells], this.Lx_d_[cells], this.x_[cells]);
            // if j == 52 {
            //     writeln("After Thomas solve for row ", j);
            //     writeln("Lx_a_: ", this.Lx_a_[cells]);
            //     writeln("Lx_b_: ", this.Lx_b_[cells]);
            //     writeln("Lx_c_: ", this.Lx_c_[cells]);
            //     writeln("Lx_d_: ", this.Lx_d_[cells]);
            //     writeln("x_ after Lx solve for row ", j, ": ", this.x_[cells]);
            // }
        }
        // writeln("Lx_a_: ", this.Lx_a_);
        // writeln("Lx_b_: ", this.Lx_b_);
        // writeln("Lx_c_: ", this.Lx_c_);
        // writeln("Lx_d_: ", this.Lx_d_);
        // writeln("x_ after Lx solve: ", this.x_);

        // Verify that A*x = b
        // var Ax : [dom_] real(64);
        // forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
        //     for i in 0..<this.spatialDisc_.mesh_.niCell_ {
        //         const cellIndex = this.spatialDisc_.mesh_.iiCell(i, j);
        //         Ax[cellIndex] = this.Lx_a_[cellIndex] * this.x_[cellIndex - 1] + this.Lx_b_[cellIndex] * this.x_[cellIndex] + this.Lx_c_[cellIndex] * this.x_[cellIndex + 1];
        //         if i == 0 && j == 52 {
        //             writeln("Checking Lx at cellIndex ", cellIndex);
        //             writeln("Lx_a_[", cellIndex, "] * x_[", cellIndex - 1, "] + Lx_b_[", cellIndex, "] * x_[", cellIndex, "] + Lx_c_[", cellIndex, "] * x_[", cellIndex + 1, "] = ",
        //             this.Lx_a_[cellIndex], " * ", this.x_[cellIndex - 1], " + ", this.Lx_b_[cellIndex], " * ", this.x_[cellIndex], " + ", this.Lx_c_[cellIndex], " * ", this.x_[cellIndex + 1], " = ", Ax[cellIndex],
        //             " (should be ", this.Lx_d_[cellIndex], ")");
        //         }
        //     }
        // }
        // writeln("Ax after Lx solve: ", Ax);
        // printArray(Ax, this.spatialDisc_.mesh_.niCell_, this.spatialDisc_.mesh_.njCell_, "Ax after Lx solve");

        this.x_ = reArrangeRow2ColumnMajor(this.x_);

        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            for i in 0..<this.spatialDisc_.mesh_.niCell_ {
                const cellIndex = this.spatialDisc_.mesh_.iiCell(i, j);

                this.Ly_d_[cellIndex] = this.x_[cellIndex];
            }
        }
        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            const cells = this.spatialDisc_.mesh_.iiCell(0, j)..this.spatialDisc_.mesh_.iiCell(this.spatialDisc_.mesh_.niCell_-1, j);
            thomasAlgorithm(this.Ly_a_[cells], this.Ly_b_[cells], this.Ly_c_[cells], this.Ly_d_[cells], this.x_[cells]);
            // if j == 52 {
            //     writeln("After Thomas solve for row ", j);
            //     writeln("Lx_a_: ", this.Lx_a_[cells]);
            //     writeln("Lx_b_: ", this.Lx_b_[cells]);
            //     writeln("Lx_c_: ", this.Lx_c_[cells]);
            //     writeln("Lx_d_: ", this.Lx_d_[cells]);
            //     writeln("x_ after Lx solve for row ", j, ": ", this.x_[cells]);
            // }
        }
        // writeln("Ly_a_: ", this.Ly_a_);
        // writeln("Ly_b_: ", this.Ly_b_);
        // writeln("Ly_c_: ", this.Ly_c_);
        // writeln("Ly_d_: ", this.Ly_d_);
        this.x_ = reArrangeRow2ColumnMajor(this.x_);
        // writeln("x_ after Ly solve: ", this.x_);

        forall j in 0..<this.spatialDisc_.mesh_.njCell_ {
            for i in 0..<this.spatialDisc_.mesh_.niCell_ {
                const cellIndex = this.spatialDisc_.mesh_.iiCell(i, j);
                const cellIndexFVM = this.spatialDisc_.meshIndex2FVMindex(i, j);

                if this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != 1 && this.spatialDisc_.cellTypesWithGhosts_[cellIndexFVM] != 0 {
                    continue;
                }
                
                // update past time level variables
                this.spatialDisc_.phiphi_m2_[cellIndexFVM] = this.spatialDisc_.phiphi_m1_[cellIndexFVM];
                
                this.spatialDisc_.rhorho_m1_[cellIndexFVM] = this.spatialDisc_.rhorho_[cellIndexFVM];
                this.spatialDisc_.beta_m1_[cellIndexFVM] = this.spatialDisc_.beta_[cellIndexFVM];
                this.spatialDisc_.phiphi_m1_[cellIndexFVM] = this.spatialDisc_.phiphi_[cellIndexFVM];
                this.spatialDisc_.uu_m1_[cellIndexFVM] = this.spatialDisc_.uu_[cellIndexFVM];
                this.spatialDisc_.vv_m1_[cellIndexFVM] = this.spatialDisc_.vv_[cellIndexFVM];

                // update solution
                this.spatialDisc_.phiphi_[cellIndexFVM] += this.inputs_.OMEGA_ * this.x_[cellIndex];
                this.spatialDisc_.phi_minus_phi_m1_[cellIndexFVM] = this.spatialDisc_.phiphi_[cellIndexFVM] - this.spatialDisc_.phiphi_m1_[cellIndexFVM];

                // if cellIndex == 27756 {
                //     writeln("Updating cell ", cellIndex, ": phi = ", this.spatialDisc_.phiphi_[cellIndexFVM],
                //     " uu = ", this.spatialDisc_.uu_[cellIndexFVM],
                //     " vv = ", this.spatialDisc_.vv_[cellIndexFVM],
                //     " rho = ", this.spatialDisc_.rhorho_[cellIndexFVM],
                //     " delta phi = ", this.x_petsc.get(cellIndex),
                //     " residual = ", this.spatialDisc_.R0_[cellIndexFVM]);
                // }

            }
        }

        return (1, "Thomas");
    }

    proc updateWakeCirculation() {
        this.spatialDisc_.computeCirculation();
        this.spatialDisc_.wakeCirculation_[0] = this.spatialDisc_.circulation_;

        forall i in this.wakeFaces_dom_ {
            const topCellFVM = this.spatialDisc_.wakeFacesTopCell_[i];
            const bottomCellFVM = this.spatialDisc_.wakeFacesBottomCell_[i];
            const u_top = this.spatialDisc_.uu_[topCellFVM];
            const u_bottom = this.spatialDisc_.uu_[bottomCellFVM];
            this.wakeVelocities_[i] = 0.5 * (u_top + u_bottom);
        }

        for i in {1..<this.wakeFaces_dom_.size} {
            const face = this.spatialDisc_.wakeFaces_[i];
            const u = this.wakeVelocities_[i];
            const dx = this.spatialDisc_.mesh_.JfacesCx_[face] - this.spatialDisc_.mesh_.JfacesCx_[face-1];
            const Gamma_i = this.spatialDisc_.wakeCirculation_[i];
            const Gamma_im1 = this.spatialDisc_.wakeCirculation_[i-1];
            this.spatialDisc_.wakeCirculation_[i] = (u / dx * (Gamma_im1) + Gamma_i / this.inputs_.CFL_) / (1.0 / this.inputs_.CFL_ + u / dx);
        }
    }

    // proc updateUnsteadyWake() {

    // }

    proc updateAngleOfAttack() {
        this.inputs_.ALPHA_ = this.inputs_.ALPHA_0_ + this.inputs_.ALPHA_AMPLITUDE_ * sin(this.inputs_.ALPHA_FREQUENCY_ * this.t_ + this.inputs_.ALPHA_PHASE_);
        this.inputs_.U_INF_ = this.inputs_.MACH_ * this.inputs_.C_INF_ * cos(this.inputs_.ALPHA_ * pi / 180.0);
        this.inputs_.V_INF_ = this.inputs_.MACH_ * this.inputs_.C_INF_ * sin(this.inputs_.ALPHA_ * pi / 180.0);
        writeln("Updated angle of attack to ", this.inputs_.ALPHA_, " degrees at time ", this.t_, " U_inf = ", this.inputs_.U_INF_, " V_inf = ", this.inputs_.V_INF_);
    }

    proc solve() {
        var time: stopwatch;
        var normalize_res0 = 1e12;

        this.spatialDisc_.initializeWakeFaces();
        this.wakeFaces_dom_ = this.spatialDisc_.wakeFaces_dom_;
        this.initializeJacobian();
        

        var its : int;
        var reason : string;
        while (normalize_res0 > this.inputs_.CONV_TOL_ && this.it_ < this.inputs_.IT_MAX_ && isNan(normalize_res0) == false) {
            this.it_ += 1;
            
            time.start();
            if this.inputs_.SOLVER_ == "gmres" {
                (its, reason) = this.eulerStep();
            } else if this.inputs_.SOLVER_ == "thomas" {
                (its, reason) = this.eulerStepThomas();
            } else {
                halt("Unknown solver type: " + this.inputs_.SOLVER_);
            }
            if this.inputs_.ALPHA_ != 0.0 {
                this.updateWakeCirculation();
            }
            const (Cl, Cd, Cm) = this.spatialDisc_.compute_aerodynamics_coefficients();

            time.stop();

            const res0 = l2Norm(this.spatialDisc_.R0_);

            if this.it_ == 1 {
                this.firstRes0 = res0;
            }

            // if this.it_ >= this.inputs_.CFL_RAMP_IT_{
            //     this.inputs_.CFL_ = min(this.inputs_.CFL_ * normalize_res0 / (res0 / firstRes0), this.inputs_.CFL_RAMP_FINAL_);
            // }

            normalize_res0 = res0 / this.firstRes0;

            writeln("t: ", time.elapsed()," Iteration: ", this.it_, " CFL: ", this.inputs_.CFL_, " Res: ", res0, 
            " Cl: ", Cl, " Cd: ", Cd, " Cm: ", Cm, " GMRES its: ", its, " Reason: ", reason, 
            " Norm Res: ", normalize_res0, " Circ: ", this.spatialDisc_.circulation_);

            this.timeList.pushBack(time.elapsed());
            this.itList.pushBack(this.it_);
            this.res0List.pushBack(res0);
            this.clList.pushBack(Cl);
            this.cdList.pushBack(Cd);
            this.cmList.pushBack(Cm);

            if this.it_ % this.inputs_.CGNS_OUTPUT_FREQ_ == 0 {
                writeln("Writing CGNS output at iteration ", this.it_);
                this.spatialDisc_.writeSolution2CGNS(timeList, itList, res0List, clList, cdList, cmList);
            }

            if its == 0 {
                // this.A_petsc.matView();
                break;
            }
        }
        this.spatialDisc_.writeSolution2CGNS(timeList, itList, res0List, clList, cdList, cmList);
    }

    proc solveUnsteady() {
        while this.t_ < this.inputs_.TIME_FINAL_ {
            this.updateAngleOfAttack();
            this.solve();
            this.t_ += this.inputs_.TIME_STEP_;
        }
    }
}

} // module fullPotentialTemporalDiscretization