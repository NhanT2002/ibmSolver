use Math;

config const MESH_FILENAME : string;
config const GEOMETRY_FILENAME : string;
config const OUTPUT_FILENAME : string;

config const X_REF : real(64);
config const Y_REF : real(64);

config const MACH : real(64);
config const ALPHA : real(64);
config const GAMMA : real(64);
config const K2 : real(64);
config const K4 : real(64);

config const CFL : real(64);
config const IT_MAX : int;
config const CONV_TOL : real(64);

record inputsConfig {
    var MESH_FILENAME_: string = MESH_FILENAME;
    var GEOMETRY_FILENAME_: string = GEOMETRY_FILENAME;
    var OUTPUT_FILENAME_: string = OUTPUT_FILENAME;

    var X_REF_ : real(64) = X_REF;
    var Y_REF_ : real(64) = Y_REF;

    var MACH_ : real(64) = MACH;
    var ALPHA_ : real(64) = ALPHA;
    var GAMMA_ : real(64) = GAMMA;
    var K2_ : real(64) = K2;
    var K4_ : real(64) = K4;

    var CFL_ : real(64) = CFL;
    var IT_MAX_ : int = IT_MAX;
    var CONV_TOL_ : real(64) = CONV_TOL;

    var RHO_INF_ : real(64) = 1.0;
    var P_INF_ : real(64) = 1.0;
    var C_INF_ : real(64) = sqrt(GAMMA_ * P_INF_ / RHO_INF_);
    var U_INF_ : real(64) = MACH_ * C_INF_* cos(ALPHA_ * (pi / 180.0));
    var V_INF_ : real(64) = MACH_ * C_INF_* sin(ALPHA_ * (pi / 180.0));
    var E_INF_ : real(64) = P_INF_ / ((GAMMA_ - 1.0) * RHO_INF_) + 0.5 * (U_INF_**2 + V_INF_**2);

    proc init() {
        writeln("MESH_FILENAME = ", MESH_FILENAME);
        writeln("GEOMETRY_FILENAME = ", GEOMETRY_FILENAME);
        writeln("OUTPUT_FILENAME = ", OUTPUT_FILENAME);

        writeln("X_REF = ", X_REF);
        writeln("Y_REF = ", Y_REF);

        writeln("MACH = ", MACH);
        writeln("ALPHA = ", ALPHA);
        writeln("GAMMA = ", GAMMA);

        writeln("K2 = ", K2);
        writeln("K4 = ", K4);

        writeln("CFL = ", CFL);
        writeln("IT_MAX = ", IT_MAX);
        writeln("CONV_TOL = ", CONV_TOL);
    }
}
