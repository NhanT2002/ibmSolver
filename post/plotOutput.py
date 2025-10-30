import h5py
import matplotlib.pyplot as plt
import numpy as np

def data(filename) :
    # Read CGNS file
    with h5py.File(filename, 'r') as f:

        X = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/x/ data'][:]
        Y = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/y/ data'][:]
        Cp = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/Cp/ data'][:]

        it = f['Base/GlobalConvergenceHistory/IterationCounters/ data'][:]
        time = f['Base/GlobalConvergenceHistory/Time/ data'][:]
        res0 = f['Base/GlobalConvergenceHistory/Res0/ data'][:]
        res1 = f['Base/GlobalConvergenceHistory/Res1/ data'][:]
        res2 = f['Base/GlobalConvergenceHistory/Res2/ data'][:]
        res3 = f['Base/GlobalConvergenceHistory/Res3/ data'][:]
        cl = f['Base/GlobalConvergenceHistory/Cl/ data'][:]
        cd = f['Base/GlobalConvergenceHistory/Cd/ data'][:]
        cm = f['Base/GlobalConvergenceHistory/Cm/ data'][:]

    return X, Y, Cp, it, time, res0, res1, res2, res3, cl, cd, cm

X, Y, Cp, it, time, res0, res1, res2, res3, cl, cd, cm = data('../output/output_test_41.cgns')

plt.figure()
plt.plot(X, Cp, "o")
plt.xlabel('x')
plt.ylabel('Cp')
plt.gca().invert_yaxis()
plt.show()

plt.figure()
plt.semilogy(it, res0, label='Res0')
plt.xlabel('Iteration')
plt.ylabel('Residuals')
plt.legend()

plt.figure()
plt.semilogy(time, res0, label='Res0')
plt.xlabel('Time (s)')
plt.ylabel('Residuals')
plt.legend()