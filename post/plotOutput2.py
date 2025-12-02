import h5py
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
import pandas as pd

def readCGNSEuler(filename) :
    # Read CGNS file
    data = {}
    with h5py.File(filename, 'r') as f:

        x = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/x/ data'][:]
        y = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/y/ data'][:]
        Cp = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/Cp/ data'][:]
        Mach = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/Mach/ data'][:]
        Density = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/Density/ data'][:]

        x_mirror = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/x_mirror/ data'][:]
        y_mirror = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/y_mirror/ data'][:]
        u = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityX/ data'][:]
        v = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityY/ data'][:]
        nx = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/nx/ data'][:]
        ny = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/ny/ data'][:]
        curvature = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/curvature/ data'][:]

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
        # res2 = f['Base/GlobalConvergenceHistory/Res2/ data'][:]
        # res3 = f['Base/GlobalConvergenceHistory/Res3/ data'][:]
        cl = f['Base/GlobalConvergenceHistory/Cl/ data'][:]
        cd = f['Base/GlobalConvergenceHistory/Cd/ data'][:]
        cm = f['Base/GlobalConvergenceHistory/Cm/ data'][:]

        Rc0 = f['Base/Zone/FLOW_SOLUTION_CC/Rc0/ data'][:]

        data['Rc0'] = Rc0

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
        # data['res2'] = res2
        # data['res3'] = res3
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
        data['curvature'] = curvature
        data['mach'] = Mach
        data['density'] = Density

    return data

def readCGNS(filename) :
    # Read CGNS file
    data = {}
    with h5py.File(filename, 'r') as f:

        x = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/x/ data'][:]
        y = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/y/ data'][:]
        Cp = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/Cp/ data'][:]
        Mach = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/Mach/ data'][:]
        Density = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/Density/ data'][:]

        x_mirror = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/x_mirror/ data'][:]
        y_mirror = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/y_mirror/ data'][:]
        u_mirror = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityX_mirror/ data'][:]
        v_mirror = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityY_mirror/ data'][:]

        u = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityX/ data'][:]
        v = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityY/ data'][:]
        nx = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/nx/ data'][:]
        ny = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/ny/ data'][:]
        curvature = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/curvature/ data'][:]

        x_ghost = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/x_ghost/ data'][:]
        y_ghost = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/y_ghost/ data'][:]
        u_ghost = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityX_ghost/ data'][:]
        v_ghost = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityY_ghost/ data'][:]
        p_ghost = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/Pressure_ghost/ data'][:]
        rho_ghost = f['WallBase/ZoneWall/WALL_FLOW_SOLUTION_NC/Density_ghost/ data'][:]

        x2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/x/ data'][:]
        y2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/y/ data'][:]
        Cp2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/Cp/ data'][:]
        Mach2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/Mach/ data'][:]
        Density2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/Density/ data'][:]

        x_mirror2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/x_mirror/ data'][:]
        y_mirror2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/y_mirror/ data'][:]
        u_mirror2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityX_mirror/ data'][:]
        v_mirror2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityY_mirror/ data'][:]

        u2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityX/ data'][:]
        v2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityY/ data'][:]
        nx2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/nx/ data'][:]
        ny2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/ny/ data'][:]

        x_ghost2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/x_ghost/ data'][:]
        y_ghost2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/y_ghost/ data'][:]
        u_ghost2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityX_ghost/ data'][:]
        v_ghost2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/VelocityY_ghost/ data'][:]
        p_ghost2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/Pressure_ghost/ data'][:]
        rho_ghost2 = f['WallBase2/ZoneWall/WALL_FLOW_SOLUTION_NC/Density_ghost/ data'][:]

        it = f['Base/GlobalConvergenceHistory/IterationCounters/ data'][:]
        time = f['Base/GlobalConvergenceHistory/Time/ data'][:]
        res0 = f['Base/GlobalConvergenceHistory/Res0/ data'][:]
        # res2 = f['Base/GlobalConvergenceHistory/Res2/ data'][:]
        # res3 = f['Base/GlobalConvergenceHistory/Res3/ data'][:]
        cl = f['Base/GlobalConvergenceHistory/Cl/ data'][:]
        cd = f['Base/GlobalConvergenceHistory/Cd/ data'][:]
        cm = f['Base/GlobalConvergenceHistory/Cm/ data'][:]

        R0 = f['Base/Zone/FLOW_SOLUTION_CC/R0/ data'][:]

        data['R0'] = R0
        data['mach_cell'] = f['Base/Zone/FLOW_SOLUTION_CC/Mach/ data'][:]

        data['x'] = x
        data['y'] = y
        data['Cp'] = Cp
        data['x_mirror'] = x_mirror
        data['y_mirror'] = y_mirror
        data['u_mirror'] = u_mirror
        data['v_mirror'] = v_mirror
        data['u'] = u
        data['v'] = v
        data['it'] = it
        data['time'] = time
        data['res0'] = res0
        data['res0'] = data['res0']/data['res0'][0]
        # data['res2'] = res2
        # data['res3'] = res3
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
        data['curvature'] = curvature
        data['mach'] = Mach
        data['density'] = Density

        data['x2'] = x2
        data['y2'] = y2
        data['Cp2'] = Cp2
        data['x_mirror2'] = x_mirror2
        data['y_mirror2'] = y_mirror2
        data['u_mirror2'] = u_mirror2
        data['v_mirror2'] = v_mirror2
        data['u2'] = u2
        data['v2'] = v2
        data['nx2'] = nx2
        data['ny2'] = ny2
        data['x_ghost2'] = x_ghost2
        data['y_ghost2'] = y_ghost2
        data['u_ghost2'] = u_ghost2
        data['v_ghost2'] = v_ghost2
        data['p_ghost2'] = p_ghost2
        data['rho_ghost2'] = rho_ghost2
        data['mach2'] = Mach2
        data['density2'] = Density2

    return data

def readGeom(filename) :
    # Read geometry file
    data = {}
    with h5py.File(filename, 'r') as f:

        x = f['Base/Geometry/GridCoordinates/CoordinateX/ data'][:]
        y = f['Base/Geometry/GridCoordinates/CoordinateY/ data'][:]
        z = f['Base/Geometry/GridCoordinates/CoordinateZ/ data'][:]

        data['x'] = x
        data['y'] = y
        data['z'] = z

    return data

data = readCGNS('../output/output_109.cgns')
data2 = readCGNSEuler('../output/output_36.cgns')
data3 = readCGNSEuler('../output/output_41.cgns')


# naca0012 = np.loadtxt('../pre/naca0012.dat')
# geom = readGeom('../pre/cylinder_geometry.cgns')
geom = readGeom('../pre/naca0012_geometry.cgns')

plt.figure()
plt.axis('equal')
plt.plot(geom['x'], geom['y'])

plt.plot(data['x'], data['y'], ".")
plt.plot(data['x_mirror'], data['y_mirror'], "x")
plt.plot(data['x_ghost'], data['y_ghost'], ".")
plt.quiver(data['x'], data['y'], data['u'], data['v'])
plt.quiver(data['x_ghost'], data['y_ghost'], data['u_ghost'], data['v_ghost'], color='g')
plt.quiver(data['x_mirror'], data['y_mirror'], data['u_mirror'], data['v_mirror'], color='orange')

# plt.quiver(data['x'], data['y'], data['nx'], data['ny'], color='r', scale=20)

plt.plot(data['x2'], data['y2'], ".")
plt.plot(data['x_mirror2'], data['y_mirror2'], ".")
plt.plot(data['x_ghost2'], data['y_ghost2'], ".")
plt.quiver(data['x2'], data['y2'], data['u2'], data['v2'])
plt.quiver(data['x_ghost2'], data['y_ghost2'], data['u_ghost2'], data['v_ghost2'], color='g')
plt.quiver(data['x_mirror2'], data['y_mirror2'], data['u_mirror2'], data['v_mirror2'], color='orange')
plt.xlim([-0.01, 0.1])
plt.show()

# for i in range(len(data['x2'])) :
#     plt.figure()
#     plt.axis('equal')
#     plt.plot(geom['x'], geom['y'])
#     # plt.plot(data['x'], data['y'], ".")
#     # plt.plot(data['x_mirror'], data['y_mirror'], "x")
#     # plt.plot(data['x_ghost'], data['y_ghost'], ".")
#     # # plt.quiver(data['x'], data['y'], data['u'], data['v'])
#     # # plt.quiver(data['x_ghost'], data['y_ghost'], data['u_ghost'], data['v_ghost'], color='g')
#     # # plt.quiver(data['x_mirror'], data['y_mirror'], 2*data['u'] - data['u_ghost'], 2*data['v'] - data['v_ghost'], color='orange')

#     # plt.quiver(data['x'], data['y'], data['nx'], data['ny'], color='r', scale=20)

#     plt.plot(data['x2'], data['y2'], ".")
#     plt.plot(data['x_mirror2'], data['y_mirror2'], ".")
#     plt.plot(data['x_ghost2'], data['y_ghost2'], ".")
#     plt.quiver(data['x2'], data['y2'], data['nx2'], data['ny2'], color='m', scale=20)
#     plt.plot(data2['Cx'], data2['Cy'], ".")
#     plt.plot(data2['Cx'][i*4:(i+1)*4], data2['Cy'][i*4:(i+1)*4], "x", markersize=10, color='red')
#     # plt.xlim([-0.5, 0.5])
#     plt.show()

plt.figure()
plt.plot(data['x'], data['Cp'], "o-")
# plt.plot(data2['x'], data2['Cp'], "s-")
plt.plot(data3['x'], data3['Cp'], "^-")
plt.xlabel('x')
plt.ylabel('Cp')
plt.gca().invert_yaxis()
plt.grid()
plt.show()

plt.figure()
plt.plot(data['x'], data['mach'], "o-")
plt.plot(data2['x'], data2['mach'], "s-")
plt.plot(data3['x'], data3['mach'], "^-")
plt.xlabel('x')
plt.ylabel('Mach')
plt.grid()
plt.show()

plt.figure()
plt.semilogy(data['it'], data['res0'], label='Res0')
plt.semilogy(data2['it'], data2['res0'], label='Res0 (data2)')
plt.semilogy(data3['it'], data3['res0'], label='Res0 (data3)')
plt.xlabel('Iteration')
plt.ylabel('Residuals')
plt.legend()
plt.grid()

plt.figure()
plt.semilogy(data['time'], data['res0'], label='Res0')
plt.semilogy(data2['time'], data2['res0'], label='Res0 (data2)')
plt.semilogy(data3['time'], data3['res0'], label='Res0 (data3)')
plt.xlabel('Time (s)')
plt.ylabel('Residual')
plt.legend()
plt.grid()

# plt.figure()
# plt.plot(data['it'], data['cl'], label='Cl')
# plt.xlabel('Iteration')
# plt.ylabel('Cl')
# plt.legend()

# plt.figure()
# plt.semilogy(data['time'], data['res0'], label='Res0')
# plt.xlabel('Time (s)')
# plt.ylabel('Residuals')
# plt.legend()
