import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('debug.txt', sep="; ")
plt.figure()
plt.axis('equal')
plt.quiver(data['x_bi'], data['y_bi'], data['u_bi'], data['v_bi'])

naca12 = np.loadtxt('pre/naca0012.dat')
plt.plot(naca12[0,:], naca12[1,:])

plt.xlim(0.2, 0.4)
plt.show()