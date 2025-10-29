import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('debug.txt', sep=", ", header=None)
x = data[0]
y = data[1]
plt.figure()
plt.plot(x, y, "o")