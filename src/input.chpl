use Math;

config const MESH_FILENAME : string;

record inputsConfig {
    var MESH_FILENAME_: string = MESH_FILENAME;

    proc init() {
        writeln("MESH_FILENAME = ", MESH_FILENAME);
    }
}
