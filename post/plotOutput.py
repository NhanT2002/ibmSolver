import h5py
import matplotlib.pyplot as plt
import numpy as np

def readCGNS(filename) :
    # Read CGNS file
    data = {}
    with h5py.File(filename, 'r') as f:

        x = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/x/ data'][:]
        y = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/y/ data'][:]
        Cp = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/Cp/ data'][:]

        x_mirror = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/x_mirror/ data'][:]
        y_mirror = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/y_mirror/ data'][:]
        u = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityX/ data'][:]
        v = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityY/ data'][:]
        nx = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/nx/ data'][:]
        ny = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/ny/ data'][:]

        x_ghost = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/x_ghost/ data'][:]
        y_ghost = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/y_ghost/ data'][:]
        u_ghost = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityX_ghost/ data'][:]
        v_ghost = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityY_ghost/ data'][:]
        p_ghost = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/Pressure_ghost/ data'][:]
        rho_ghost = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/Density_ghost/ data'][:]

        it = f['Base/GlobalConvergenceHistory/IterationCounters/ data'][:]
        time = f['Base/GlobalConvergenceHistory/Time/ data'][:]
        res0 = f['Base/GlobalConvergenceHistory/Res0/ data'][:]
        res1 = f['Base/GlobalConvergenceHistory/Res1/ data'][:]
        res2 = f['Base/GlobalConvergenceHistory/Res2/ data'][:]
        res3 = f['Base/GlobalConvergenceHistory/Res3/ data'][:]
        cl = f['Base/GlobalConvergenceHistory/Cl/ data'][:]
        cd = f['Base/GlobalConvergenceHistory/Cd/ data'][:]
        cm = f['Base/GlobalConvergenceHistory/Cm/ data'][:]

        data['x'] = x
        data['y'] = y
        data['Cp'] = Cp
        data['x_mirror'] = x_mirror
        data['y_mirror'] = y_mirror
        data['u'] = u
        data['v'] = v
        data['it'] = it
        data['time'] = time
        data['res0'] = res0
        data['res1'] = res1
        data['res2'] = res2
        data['res3'] = res3
        data['cl'] = cl
        data['cd'] = cd
        data['cm'] = cm
        data['nx'] = nx
        data['ny'] = ny
        data['x_ghost'] = x_ghost
        data['y_ghost'] = y_ghost
        data['u_ghost'] = u_ghost
        data['v_ghost'] = v_ghost
        data['p_ghost'] = p_ghost
        data['rho_ghost'] = rho_ghost

    return data

data = readCGNS('../output/output_24.cgns')

naca0012 = np.loadtxt('naca0012.dat')

plt.figure()
plt.axis('equal')
plt.plot(naca0012[0, :], naca0012[1, :])
plt.plot(data['x'], data['y'], ".")
plt.plot(data['x_mirror'], data['y_mirror'], "x")
plt.plot(data['x_ghost'], data['y_ghost'], ".")
plt.quiver(data['x'], data['y'], data['u'], data['v'])
# plt.quiver(data['x'], data['y'], data['nx'], data['ny'], color='r', scale=20)
plt.quiver(data['x_ghost'], data['y_ghost'], data['u_ghost'], data['v_ghost'], color='g')
plt.quiver(data['x_mirror'], data['y_mirror'], 2*data['u'] - data['u_ghost'], 2*data['v'] - data['v_ghost'], color='orange')
plt.xlim([0.0, 0.1])

plt.figure()
plt.plot(data['x'], data['Cp'], "o")
plt.xlabel('x')
plt.ylabel('Cp')
plt.gca().invert_yaxis()
plt.show()

plt.figure()
plt.semilogy(data['it'], data['res0'], label='Res0')
plt.xlabel('Iteration')
plt.ylabel('Residuals')
plt.legend()

# plt.figure()
# plt.semilogy(data['time'], data['res0'], label='Res0')
# plt.xlabel('Time (s)')
# plt.ylabel('Residuals')
# plt.legend()
