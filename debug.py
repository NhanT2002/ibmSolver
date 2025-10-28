import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('debug.txt', sep=", ")
plt.figure()
plt.axis('equal')
plt.plot(data['xCells'], data['yCells'], 'o')
# plt.quiver(data['xCells'], data['yCells'], data['nx_bi'], data['ny_bi'])
plt.plot(data['x_mirror'], data['y_mirror'],'x', color='red')
plt.plot(data['x_bi'], data['y_bi'],'x', color='orange')

naca12 = np.loadtxt('pre/naca0012.dat')
plt.plot(naca12[0,:], naca12[1,:])

plt.xlim(0.2, 0.4)