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

config const MACH : real(64);
config const ALPHA : real(64);
config const GAMMA : real(64);
config const K2 : real(64);
config const K4 : real(64);

config const CFL : real(64);
config const CFL_RAMP_FACTOR : real(64);
config const CFL_RAMP_IT : int;
config const CFL_RAMP_MAX : real(64);
config const OMEGA : real(64);
config const IT_MAX : int;
config const CONV_TOL : real(64);
config const RESIDUAL_SMOOTHING : bool;

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

    var MACH_ : real(64) = MACH;
    var ALPHA_ : real(64) = ALPHA;
    var GAMMA_ : real(64) = GAMMA;
    var K2_ : real(64) = K2;
    var K4_ : real(64) = K4;

    var CFL_ : real(64) = CFL;
    var CFL_RAMP_FACTOR_ : real(64) = CFL_RAMP_FACTOR;
    var CFL_RAMP_IT_ : int = CFL_RAMP_IT;
    var CFL_RAMP_MAX_ : real(64) = CFL_RAMP_MAX;
    var OMEGA_ : real(64) = OMEGA;
    var IT_MAX_ : int = IT_MAX;
    var CONV_TOL_ : real(64) = CONV_TOL;
    var RESIDUAL_SMOOTHING_ : bool = RESIDUAL_SMOOTHING;

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

        writeln("K2 = ", K2);
        writeln("K4 = ", K4);

        writeln("CFL = ", CFL);
        writeln("CFL_RAMP_FACTOR = ", CFL_RAMP_FACTOR);
        writeln("CFL_RAMP_IT = ", CFL_RAMP_IT);
        writeln("CFL_RAMP_MAX = ", CFL_RAMP_MAX);
        writeln("OMEGA = ", OMEGA);
        writeln("IT_MAX = ", IT_MAX);
        writeln("CONV_TOL = ", CONV_TOL);
        writeln("RESIDUAL_SMOOTHING = ", RESIDUAL_SMOOTHING);

        writeln("GMRES_RTOL = ", GMRES_RTOL);
        writeln("GMRES_ATOL = ", GMRES_ATOL);
        writeln("GMRES_DTOL = ", GMRES_DTOL);
        writeln("GMRES_MAXIT = ", GMRES_MAXIT);
        writeln("GMRES_RESTART = ", GMRES_RESTART);



        if FLOW == "euler" {
            RHO_INF_ = 1.0;
            P_INF_ = 1.0;
            C_INF_ = sqrt(GAMMA_ * P_INF_ / RHO_INF_);
            U_INF_ = MACH_ * C_INF_* cos(ALPHA_ * (pi / 180.0));
            V_INF_ = MACH_ * C_INF_* sin(ALPHA_ * (pi / 180.0));
            E_INF_ = P_INF_ / ((GAMMA_ - 1.0) * RHO_INF_) + 0.5 * (U_INF_**2 + V_INF_**2);
            Q_INF_ = 0.5 * RHO_INF_ * (U_INF_**2 + V_INF_**2);
            S_REF_ = 1.0; // Reference area
            C_REF_ = 1.0; // Reference chord
        }
        else {
            RHO_INF_ = 1.0;
            P_INF_ = RHO_INF_**GAMMA_ / (GAMMA_ * MACH_**2);
            C_INF_ = sqrt((RHO_INF_**(GAMMA_ - 1.0)) / MACH_**2);
            U_INF_ = MACH_ * C_INF_* cos(ALPHA_ * (pi / 180.0));
            V_INF_ = MACH_ * C_INF_* sin(ALPHA_ * (pi / 180.0));
            Q_INF_ = 0.5 * RHO_INF_ * (U_INF_**2 + V_INF_**2);
            S_REF_ = 1.0; // Reference area
            C_REF_ = 1.0; // Reference chord
        }
    }
}
