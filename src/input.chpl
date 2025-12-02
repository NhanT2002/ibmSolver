use Math;
config const FLOW : string;

config const MESH_FILENAME : string;
config const GEOMETRY_FILENAME : string;
config const OUTPUT_FILENAME : string;
config const CGNS_OUTPUT_FREQ : int;

config const INITIAL_SOLUTION : string = "";

config const X_REF : real(64);
config const Y_REF : real(64);

config const X_LE : real(64);
config const Y_LE : real(64);
config const X_TE : real(64);
config const Y_TE : real(64);

config const N_KLS : int;

config const MACH : real(64);
config const ALPHA : real(64);
config const GAMMA : real(64);
config const K2 : real(64);
config const K4 : real(64);

config const ALPHA_0 : real(64);
config const ALPHA_AMPLITUDE : real(64);
config const ALPHA_FREQUENCY : real(64);
config const ALPHA_PHASE : real(64);
config const TIME_STEP : real(64);
config const TIME_FINAL : real(64);

config const CFL : real(64);
config const CFL_RAMP_FACTOR : real(64);
config const CFL_RAMP_IT : int;
config const CFL_RAMP_FINAL : real(64);
config const OMEGA : real(64);
config const IT_MAX : int;
config const CONV_TOL : real(64);
config const RESIDUAL_SMOOTHING : bool;

config const SOLVER : string;

config const GMRES_PRECON : string;
config const GMRES_RTOL : real(64);
config const GMRES_ATOL : real(64);
config const GMRES_DTOL : real(64);
config const GMRES_MAXIT : int;
config const GMRES_RESTART : int;

record inputsConfig {
    var FLOW_: string = FLOW;

    var MESH_FILENAME_: string = MESH_FILENAME;
    var GEOMETRY_FILENAME_: string = GEOMETRY_FILENAME;
    var OUTPUT_FILENAME_: string = OUTPUT_FILENAME;
    var CGNS_OUTPUT_FREQ_: int = CGNS_OUTPUT_FREQ;

    var INITIAL_SOLUTION_ : string = INITIAL_SOLUTION;

    var X_REF_ : real(64) = X_REF;
    var Y_REF_ : real(64) = Y_REF;

    var X_LE_ : real(64) = X_LE;
    var Y_LE_ : real(64) = Y_LE;
    var X_TE_ : real(64) = X_TE;
    var Y_TE_ : real(64) = Y_TE;

    var N_KLS_ : int = N_KLS; // Number of cells in k-exact least squares stencil

    var MACH_ : real(64) = MACH;
    var MACH_2_ : real(64) = MACH * MACH;
    var ALPHA_ : real(64) = ALPHA;
    var GAMMA_ : real(64) = GAMMA;
    var K2_ : real(64) = K2;
    var K4_ : real(64) = K4;

    var ALPHA_0_: real(64) = ALPHA_0;
    var ALPHA_AMPLITUDE_: real(64) = ALPHA_AMPLITUDE;
    var ALPHA_FREQUENCY_: real(64) = ALPHA_FREQUENCY;
    var ALPHA_PHASE_: real(64) = ALPHA_PHASE;
    var TIME_STEP_ : real(64) = TIME_STEP;
    var TIME_FINAL_ : real(64) = TIME_FINAL;

    var CFL_ : real(64) = CFL;
    var CFL_RAMP_FACTOR_ : real(64) = CFL_RAMP_FACTOR;
    var CFL_RAMP_IT_ : int = CFL_RAMP_IT;
    var CFL_RAMP_FINAL_ : real(64) = CFL_RAMP_FINAL;
    var OMEGA_ : real(64) = OMEGA;
    var IT_MAX_ : int = IT_MAX;
    var CONV_TOL_ : real(64) = CONV_TOL;
    var RESIDUAL_SMOOTHING_ : bool = RESIDUAL_SMOOTHING;

    var SOLVER_ : string = SOLVER;

    var GMRES_PRECON_ : string = GMRES_PRECON;
    var GMRES_RTOL_: real(64) = GMRES_RTOL;
    var GMRES_ATOL_: real(64) = GMRES_ATOL;
    var GMRES_DTOL_: real(64) = GMRES_DTOL;
    var GMRES_MAXIT_: int = GMRES_MAXIT;
    var GMRES_RESTART_: int = GMRES_RESTART;    

    
    var RHO_INF_ : real(64) = 1.0;
    var P_INF_ : real(64) = 1.0;
    var C_INF_ : real(64) = sqrt(GAMMA_ * P_INF_ / RHO_INF_);
    var U_INF_ : real(64) = MACH_ * C_INF_* cos(ALPHA_ * (pi / 180.0));
    var V_INF_ : real(64) = MACH_ * C_INF_* sin(ALPHA_ * (pi / 180.0));
    var E_INF_ : real(64) = P_INF_ / ((GAMMA_ - 1.0) * RHO_INF_) + 0.5 * (U_INF_**2 + V_INF_**2);
    var Q_INF_ : real(64) = 0.5 * RHO_INF_ * (U_INF_**2 + V_INF_**2);
    var S_REF_ : real(64) = 1.0; // Reference area
    var C_REF_ : real(64) = 1.0; // Reference chord

    proc init() {
        writeln("----- Input Configuration -----");
        writeln("FLOW = ", FLOW);
        writeln("MESH_FILENAME = ", MESH_FILENAME);
        writeln("GEOMETRY_FILENAME = ", GEOMETRY_FILENAME);
        writeln("OUTPUT_FILENAME = ", OUTPUT_FILENAME);
        writeln("CGNS_OUTPUT_FREQ = ", CGNS_OUTPUT_FREQ);

        writeln("INITIAL_SOLUTION = ", INITIAL_SOLUTION);

        writeln("X_REF = ", X_REF);
        writeln("Y_REF = ", Y_REF);

        writeln("MACH = ", MACH);
        writeln("ALPHA = ", ALPHA);
        writeln("GAMMA = ", GAMMA);

        writeln("ALPHA_0 = ", ALPHA_0);
        writeln("ALPHA_AMPLITUDE = ", ALPHA_AMPLITUDE);
        writeln("ALPHA_FREQUENCY = ", ALPHA_FREQUENCY);
        writeln("ALPHA_PHASE = ", ALPHA_PHASE);

        writeln("K2 = ", K2);
        writeln("K4 = ", K4);

        writeln("CFL = ", CFL);
        writeln("CFL_RAMP_FACTOR = ", CFL_RAMP_FACTOR);
        writeln("CFL_RAMP_IT = ", CFL_RAMP_IT);
        writeln("CFL_RAMP_FINAL = ", CFL_RAMP_FINAL);
        writeln("OMEGA = ", OMEGA);
        writeln("IT_MAX = ", IT_MAX);
        writeln("CONV_TOL = ", CONV_TOL);
        writeln("RESIDUAL_SMOOTHING = ", RESIDUAL_SMOOTHING);

        writeln("SOLVER = ", SOLVER);

        writeln("GMRES_RTOL = ", GMRES_RTOL);
        writeln("GMRES_ATOL = ", GMRES_ATOL);
        writeln("GMRES_DTOL = ", GMRES_DTOL);
        writeln("GMRES_MAXIT = ", GMRES_MAXIT);
        writeln("GMRES_RESTART = ", GMRES_RESTART);
    }
}

proc ref inputsConfig.initializeFlowField() {
    // Initial flow field based on freestream conditions
    if this.FLOW_ == "euler" {
        this.RHO_INF_ = 1.0;
        this.P_INF_ = 1.0;
        this.C_INF_ = sqrt(this.GAMMA_ * this.P_INF_ / this.RHO_INF_);
        this.U_INF_ = this.MACH_ * this.C_INF_* cos(this.ALPHA_ * (pi / 180.0));
        this.V_INF_ = this.MACH_ * this.C_INF_* sin(this.ALPHA_ * (pi / 180.0));
        this.E_INF_ = this.P_INF_ / ((this.GAMMA_ - 1.0) * this.RHO_INF_) + 0.5 * (this.U_INF_**2 + this.V_INF_**2);
        this.Q_INF_ = 0.5 * this.RHO_INF_ * (this.U_INF_**2 + this.V_INF_**2);
        this.S_REF_ = 1.0; // Reference area
        this.C_REF_ = 1.0; // Reference chord
    }
    else if this.FLOW_ == "fullPotential" {
        this.RHO_INF_ = 1.0;
        this.C_INF_ = 1.0;
        this.P_INF_ = this.RHO_INF_ * this.C_INF_**2 / this.GAMMA_;
        this.U_INF_ = this.MACH_ * this.C_INF_* cos(this.ALPHA_ * (pi / 180.0));
        this.V_INF_ = this.MACH_ * this.C_INF_* sin(this.ALPHA_ * (pi / 180.0));
        this.Q_INF_ = 0.5 * this.RHO_INF_ * (this.U_INF_**2 + this.V_INF_**2);
        this.S_REF_ = 1.0; // Reference area
        this.C_REF_ = 1.0; // Reference chord
    }
    else {
        this.RHO_INF_ = 1.0;
        this.C_INF_ = 1.0 / this.MACH_;
        this.P_INF_ = this.RHO_INF_**this.GAMMA_ / (this.GAMMA_ * this.MACH_**2);
        this.U_INF_ = this.MACH_ * this.C_INF_* cos(this.ALPHA_ * (pi / 180.0));
        this.V_INF_ = this.MACH_ * this.C_INF_* sin(this.ALPHA_ * (pi / 180.0));
        this.Q_INF_ = 0.5 * this.RHO_INF_ * (this.U_INF_**2 + this.V_INF_**2);
        this.S_REF_ = 1.0; // Reference area
        this.C_REF_ = 1.0;
    }
}
