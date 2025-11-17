import h5py
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

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

data = readCGNS('../output/output_198.cgns')
data2 = readCGNS('../output/output_239.cgns')

# data = readCGNS('../output/mesh_0-01/M05_A125.cgns')

# naca0012 = np.loadtxt('../pre/naca0012.dat')
geom = readGeom('../pre/cylinder_geometry.cgns')

# Calculate curvature
# x = sp.symbols('x')
# t = 0.12
# f = 5 * t * (0.2969 * sp.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)
# f_prime = sp.diff(f, x)
# f_double_prime = sp.diff(f_prime, x)
# curvature = sp.simplify(sp.Abs(f_double_prime) / (1 + f_prime**2)**(3/2))
# print("Curvature expression:")
# sp.pprint(curvature)

# # Evaluate curvature at discrete points
# x_vals = np.linspace(0, 1, 1001)
# curvature_func = sp.lambdify(x, curvature, 'numpy')
# curvature_vals = curvature_func(x_vals)

# plt.figure()
# plt.semilogy(x_vals, curvature_vals, label='Curvature of NACA 0012')
# plt.semilogy(data['x'], data['curvature'], 'ro', markersize=2, label='Computed Curvature from CGNS')
# plt.xlabel('x')
# plt.ylabel('Curvature')
# plt.title('Curvature Distribution along NACA 0012 Airfoil')
# plt.legend()
# plt.grid()
# plt.show()

plt.figure()
plt.axis('equal')
plt.plot(geom['x'], geom['y'])
plt.plot(data['x'], data['y'], ".")
plt.plot(data['x_mirror'], data['y_mirror'], "x")
plt.plot(data['x_ghost'], data['y_ghost'], ".")
plt.quiver(data['x'], data['y'], data['u'], data['v'])
# plt.quiver(data['x'], data['y'], data['nx'], data['ny'], color='r', scale=20)
plt.quiver(data['x_ghost'], data['y_ghost'], data['u_ghost'], data['v_ghost'], color='g')
plt.quiver(data['x_mirror'], data['y_mirror'], 2*data['u'] - data['u_ghost'], 2*data['v'] - data['v_ghost'], color='orange')
# plt.xlim([-0.5, 0.5])

plt.figure()
plt.plot(data['x'], data['Cp'], "o-")
plt.plot(data2['x'], data2['Cp'], "x-")
plt.xlabel('x')
plt.ylabel('Cp')
plt.gca().invert_yaxis()
plt.grid()
plt.show()

plt.figure()
plt.plot(data['x'], data['mach'], "o-")
plt.plot(data2['x'], data2['mach'], "x-")
plt.xlabel('x')
plt.ylabel('Mach')
plt.grid()
plt.show()

plt.figure()
plt.semilogy(data['it'], data['res0'], label='Res0')
plt.semilogy(data2['it'], data2['res0'], label='Res0 (later)')
plt.xlabel('Iteration')
plt.ylabel('Residuals')
plt.legend()
plt.grid()

plt.figure()
plt.semilogy(data['time'], data['res0'], label='Res0')
plt.semilogy(data2['time'], data2['res0'], label='Res0 (later)')
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
