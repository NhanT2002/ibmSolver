module eikonal {

use mesh;
use input;inputsConfig;
use Math;
use IO;
use Heap;
use List;
use Set;
use Map;
use Sort;

class eikonal {
    var mesh_ : shared meshData;
    var inputs_ : inputsConfig;

    var narrowBand_dom_ : domain(1) = {1..0};
    var narrowBand_ : [narrowBand_dom_] int;

    proc init(ref mesh : shared meshData, ref inputs : inputsConfig) {
        this.mesh_ = mesh;
        this.inputs_ = inputs;
    }

    proc init(ref mesh : shared meshData, ref inputs : inputsConfig, narrowBand : [] int) {
        this.mesh_ = mesh;
        this.inputs_ = inputs;
        this.narrowBand_dom_ = {0..narrowBand.size-1};
        this.narrowBand_ = narrowBand;
    }

    proc initalizeNarrowBand() {
        var narrowBandSet = new set(int);
        // Identify narrow band cells (fluid cells adjacent to solid cells)
        for j in 0..<this.mesh_.njCell_ {
            for i in 0..<this.mesh_.niCell_ {
                var idx = this.mesh_.iiCell(i, j);
                if this.mesh_.cellTypes_[idx] == 0 {
                    // solid cell, check neighbors
                    var neighbors = [idx - 1, idx + 1, idx - this.mesh_.niCell_, idx + this.mesh_.niCell_];
                    for nidx in neighbors {
                        if this.mesh_.cellTypes_[nidx] == 1 {
                            // neighbor is fluid cell, add to narrow band
                            narrowBandSet.add(nidx);
                        }
                    }
                }
            }
        }
        this.narrowBand_dom_ = {0..narrowBandSet.size-1};
        this.narrowBand_ = narrowBandSet.toArray();
        // writeln("Narrow band initialized with ", this.narrowBand_.size, " cells.");
        // writeln("Narrow band cells: ", this.narrowBand_);

        // var narrowBandValues : [this.narrowBand_dom_] real;
        // var narrowBandValues2Idx : map(real, int);
        // // Initialize distance values in narrow band
        // for idx in this.narrowBand_dom_ {
        //     var cellIdx = this.narrowBand_[idx];
        //     const x = this.mesh_.xCells_[cellIdx];
        //     const y = this.mesh_.yCells_[cellIdx];
        //     narrowBandValues[idx] = x**2 + y**2; // Example: distance from origin
        //     narrowBandValues2Idx[narrowBandValues[idx]] = cellIdx;
        // }

        // var dist = this.mesh_.phiCells_;

        // record comp : relativeComparator { }

        // proc comp.compare(x, y) {
        //     const x_idx = x: int;
        //     const y_idx = y: int;
        //     const x_value = dist[x_idx];
        //     const y_value = dist[y_idx];
        //     return y_value - x_value;
        // }

        // var reverseComparator: comp;

        // var minHeap = createHeap(this.narrowBand_, false, comparator = reverseComparator);
        // writeln("Heap created with ", minHeap.size, " elements.");
        // writeln("minHeap = ", minHeap);
        // this.printHeap(minHeap);

        // writeln("minHeap.top() = ", minHeap.top(), " value = ", dist[minHeap.top()]);
        // minHeap.pop();
        // writeln("After pop, minHeap.top() = ", minHeap.top(), " value = ", dist[minHeap.top()]);
        // this.printHeap(minHeap);

    }

    proc solve2(ref T : [] real, const ref xCells : [] real, const ref yCells : [] real, F2 : real, cell : int) {
        const leftCell = cell - 1;
        const rightCell = cell + 1;
        const bottomCell = cell - this.mesh_.niCellWithGhosts_;
        const topCell = cell + this.mesh_.niCellWithGhosts_;

        const T1 = min(T[leftCell], T[rightCell]);
        const T2 = min(T[bottomCell], T[topCell]);

        var dx = 0.0;
        var dy = 0.0;
        if T1 == T[leftCell] {
            dx = xCells[cell] - xCells[leftCell];
        }
        else {
            dx = xCells[rightCell] - xCells[cell];
        }
        if T2 == T[bottomCell] {
            dy = yCells[cell] - yCells[bottomCell];
        }
        else {
            dy = yCells[topCell] - yCells[cell];
        }
        const a = 1.0/dx**2 + 1.0/dy**2;
        const b = -2.0*(T1/dx**2 + T2/dy**2);
        const c = T1**2/dx**2 + T2**2/dy**2 - F2;

        const T_new = this.solveQuadratic(a, b, c);
        // check if T_new is greater than max(T1, T2)
        if T_new >= max(T1, T2) {
            T[cell] = T_new;
        }
        else {
            T[cell] = min(T1 + abs(dx)*sqrt(F2), T2 + abs(dy)*sqrt(F2));
        }

        // if cell == 5865 {
        //     writeln("Eikonal solve at cell ", cell, ":", " F2 = ", F2);
        //     writeln(" T1 = ", T1, ", T2 = ", T2);
        //     writeln(" a = ", a, ", b = ", b, ", c = ", c);
        //     writeln(" T_new = ", T_new, ", updated T[cell] = ", T[cell]);
        //     writeln("|grad T|^2 = ", ((T[cell] - T1)/dx)**2 + ((T[cell] - T2)/dy)**2);
        // }
    }

    proc solve(ref T : [] real, const ref xCells : [] real, const ref yCells : [] real, F2 : real, cell : int) {
        // Solve the Eikonal equation at cellIdx using upwind scheme
        const leftCell = cell - 1;
        const leftCellm1 = cell - 2;
        const rightCell = cell + 1;
        const rightCellp1 = cell + 2;
        const bottomCell = cell - this.mesh_.niCellWithGhosts_;
        const bottomCellm1 = cell - 2 * this.mesh_.niCellWithGhosts_;
        const topCell = cell + this.mesh_.niCellWithGhosts_;
        const topCellp1 = cell + 2 * this.mesh_.niCellWithGhosts_;

        // backward differences
        const x0_back = xCells[leftCellm1];
        const x1_back = xCells[leftCell];
        const x2_back = xCells[cell];
        const x0_coeff_back = (x2_back - x1_back)/((x0_back - x1_back)*(x0_back - x2_back));
        const x1_coeff_back = (x2_back - x0_back)/((x1_back - x0_back)*(x1_back - x2_back));
        const x2_coeff_back = (2*x2_back - x0_back - x1_back)/((x2_back - x0_back)*(x2_back - x1_back));
        const T1_back = (T[leftCellm1]*x0_coeff_back + T[leftCell]*x1_coeff_back) / x2_coeff_back;
        // forward differences
        const x0_forw = xCells[cell];
        const x1_forw = xCells[rightCell];
        const x2_forw = xCells[rightCellp1];
        const x0_coeff_forw = (2*x0_forw - x1_forw - x2_forw)/((x0_forw - x1_forw)*(x0_forw - x2_forw));
        const x1_coeff_forw = (x0_forw - x2_forw)/((x1_forw - x0_forw)*(x1_forw - x2_forw));
        const x2_coeff_forw = (x0_forw - x1_forw)/((x2_forw - x0_forw)*(x2_forw - x1_forw));
        const T1_forw = (T[rightCell]*x1_coeff_forw + T[rightCellp1]*x2_coeff_forw) / x0_coeff_forw;

        const T1 = max(T1_back, T1_forw, 0);

        // backward differences
        const y0_back = yCells[bottomCellm1];
        const y1_back = yCells[bottomCell];
        const y2_back = yCells[cell];
        const y0_coeff_back = (y2_back - y1_back)/((y0_back - y1_back)*(y0_back - y2_back));
        const y1_coeff_back = (y2_back - y0_back)/((y1_back - y0_back)*(y1_back - y2_back));
        const y2_coeff_back = (2*y2_back - y0_back - y1_back)/((y2_back - y0_back)*(y2_back - y1_back));
        const T2_back = (T[bottomCellm1]*y0_coeff_back + T[bottomCell]*y1_coeff_back) / y2_coeff_back;
        // forward differences
        const y0_forw = yCells[cell];
        const y1_forw = yCells[topCell];
        const y2_forw = yCells[topCellp1];
        const y0_coeff_forw = (2*y0_forw - y1_forw - y2_forw)/((y0_forw - y1_forw)*(y0_forw - y2_forw));
        const y1_coeff_forw = (y0_forw - y2_forw)/((y1_forw - y0_forw)*(y1_forw - y2_forw));
        const y2_coeff_forw = (y0_forw - y1_forw)/((y2_forw - y0_forw)*(y2_forw - y1_forw));
        const T2_forw = (T[topCell]*y1_coeff_forw + T[topCellp1]*y2_coeff_forw) / y0_coeff_forw;

        const T2 = max(T2_back, T2_forw, 0);

        
        var a = 0.0;
        var b = 0.0;
        var c = 0.0;
        var x1 = 0.0;
        var y2 = 0.0;
        if T1 == T1_back {
            a += x2_coeff_back**2;
            b += 2*x2_coeff_back*T1;
            c += T1**2;
            x1 = x1_back;
        }
        else {
            a += x0_coeff_forw**2;
            b += 2*x0_coeff_forw*T1;
            c += T1**2;
            x1 = x1_forw;
        }
        if T2 == T2_back {
            a += y2_coeff_back**2;
            b += 2*y2_coeff_back*T2;
            c += T2**2;
            y2 = y1_back;
        }
        else {
            a += y0_coeff_forw**2;
            b += 2*y0_coeff_forw*T2;
            c += T2**2;
            y2 = y1_forw;
        }
        c -= F2;

        const T_new = this.solveQuadratic(a, b, c);
        // check if T_new is greater than max(T1, T2)
        if T_new >= max(T1, T2) {
            T[cell] = T_new;
        }
        else {
            T[cell] = min(T1 + abs(xCells[cell]-x1)*sqrt(F2), T2 + abs(yCells[cell]-y2)*sqrt(F2));
        }
        if cell == 6084 {
            writeln("Eikonal solve at cell ", cell, ":", " F2 = ", F2);
            writeln(" T1 = ", T1, ", T2 = ", T2);
            writeln(" a = ", a, ", b = ", b, ", c = ", c);
            writeln(" T_new = ", T_new, ", updated T[cell] = ", T[cell]);
        }
    }

    proc solveQuadratic(a : real, b : real, c : real) : real {
        const disc = b**2 - 4.0 * a * c;
        if disc < 0.0 {
            return -1.0; // No real solution
        }
        const sqrt_disc = sqrt(disc);
        const sol1 = (-b + sqrt_disc) / (2.0 * a);
        const sol2 = (-b - sqrt_disc) / (2.0 * a);
        
        return max(sol1, sol2);
    }

    proc printHeap(minHeap: heap(int), const ref values: [] real) {
        // print root and its children
        const heap_array = minHeap.toArray();
        const n = minHeap.size;
        var level = 0;
        var start_index = 0;
        while start_index < n {
            const elements_on_level = 2**level;
            const end_index = min(start_index + elements_on_level, n);
            const current_level_elements = heap_array[start_index..end_index-1];
            writeln("Level ", level, ": ", current_level_elements, " values = ", values[current_level_elements]);
            start_index = end_index;
            level += 1;
        }
    }
}

} // module eikonal