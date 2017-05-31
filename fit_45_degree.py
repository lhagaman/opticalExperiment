from fitting_library import plot_data
import numpy as np

data = np.loadtxt("opticalExperiment/45_degree_IBA_in_air.txt", skiprows=1)

xdata = [d[0] for d in data]
ydata = [d[1] for d in data]

# if the lines from back to front annoy you:
xdata_partial = xdata[0:len(xdata) / 5]
ydata_partial = ydata[0:len(ydata) / 5]

# air
n_0 = 1
# manually chosen, would normally be known during experiment,
# but I changed it because error on this angle was large
theta_i = 54
# phi_r=0 for all measurements we'll make
phi_r = 0

plot_data(xdata, ydata, n_0, theta_i, phi_r, "45 degrees IBA in air")
