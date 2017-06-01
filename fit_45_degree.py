from fitting_library import plot_data_with_fit
import numpy as np

data = np.loadtxt("opticalExperiment/45_degree_IBA_in_air.txt", skiprows=1)

theta_r_array = [d[0] for d in data]
intensity_array = [d[1] for d in data]

"""
# if the lines from back to front annoy you:
theta_r_array = theta_r_array[0:len(theta_r_array) / 5]
intensity_array = intensity_array[0:len(intensity_array) / 5]
"""

numpoints = len(theta_r_array)

# air
n_0 = 1

# manually chosen, would normally be known during experiment,
# but I changed it because error on this angle was large
theta_i = 54

# phi_r=0 for all measurements we'll make
phi_r = 0

points = []
for i in range(numpoints):
    theta_r = theta_r_array[i]
    intensity = intensity_array[i]
    points.append([theta_r, phi_r, theta_i, n_0, intensity])

plot_data_with_fit(points, "45 degrees IBA in air")
