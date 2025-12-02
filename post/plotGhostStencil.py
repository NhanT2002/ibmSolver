import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import random

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

data = pd.read_csv('../output/output.cgns_ghostCellNeighbors.txt')

geom = readGeom('../pre/naca0012_geometry.cgns')

plt.figure()
plt.axis('equal')
plt.plot(geom['x'], geom['y'])
colors = plt.cm.hsv(np.linspace(0, 1, len(data)))
for i in [0, len(data)-1] :
    color = random.choice(colors)
    cx = data.iloc[i]['ghostCellCx']
    cy = data.iloc[i]['ghostCellCy']
    neighborsCx = np.array(data.iloc[i]['neighborCellCx'].split(' '), dtype=float)
    neighborsCy = np.array(data.iloc[i]['neighborCellCy'].split(' '), dtype=float)

    plt.plot(cx, cy, 'o', color=color)
    plt.plot(neighborsCx, neighborsCy, 'o', color=color, markersize=4)
plt.xlim([0.8, 1.1])