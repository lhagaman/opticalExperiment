from fitting_library import plot_data
import numpy as np

data = np.loadtxt("45_degree_IBA_in_air.txt", skiprows=1)

xdata = [d[0] for d in data]
ydata = [d[1] for d in data]
xdata_partial = xdata[0:len(xdata) / 5]
ydata_partial = ydata[0:len(ydata) / 5]

n_0 = 1
theta_i = 54
phi_r = 0

plot_data(xdata, ydata, theta_i, n_0, "45 degrees IBA in air")
