import numpy as np
import matplotlib.pyplot as plt

from TSTR_fit import BRIDF_plotter

x = np.linspace(0, 90, 1000)
y_45 = BRIDF_plotter(x, 0., 45., 1.5, 0.5, [0, 1.35, 0.3])
y_60 = BRIDF_plotter(x, 0., 60., 1.5, 0.5, [0, 1.35, 0.3])
y_75 = BRIDF_plotter(x, 0., 75., 1.5, 0.5, [0, 1.35, 0.3])

plt.plot(x, y_45, label="45 degrees")
plt.plot(x, y_60, label="60 degrees")
plt.plot(x, y_75, label="75 degrees")

plt.legend()
plt.title("n_0=1.5, rho_l=0, gamma=0.3, n=1.35")

plt.show()
