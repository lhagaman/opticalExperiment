import torrance_sparrow_fit
import numpy as np
import matplotlib.pyplot as plt

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
theta_i_in_degrees = 54

# phi_r=0 for all measurements we'll make
phi_r_in_degrees = 0

polarization = 0

for polarization in [0, 1, 2]:
    plt.figure(polarization)
    if (polarization == 0):
        p = "unpolarized"
    elif (polarization == 1):
        p = "s-polarized"
    elif (polarization == 2):
        p = "p-polarized"
    title = "45 degree IBA in air, Torrance-Sparrow Fit, " + p

    # each point has form [theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity]
    points = []
    for i in range(numpoints):
        theta_r = theta_r_array[i]
        intensity = intensity_array[i]
        points.append([theta_r, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity])

    theta_r_in_degrees_array = []
    intensity_array = []
    for point in points:
        theta_r_in_degrees_array.append(point[0])
        intensity_array.append(point[5])

    parameters = torrance_sparrow_fit.fit_parameters(points)

    rho_L = parameters[0]
    n = parameters[1]
    gamma = parameters[2]

    fitted_data = torrance_sparrow_fit.BRIDF_plotter(
        theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, rho_L, n, gamma, polarization)

    plt.plot(theta_r_in_degrees_array, intensity_array, label="experimental")
    # plt.plot(theta_r_in_degrees_array, [x[0] for x in fitted_data], label="fitted specular")
    # plt.plot(theta_r_in_degrees_array, [x[1] for x in fitted_data], label="fitted diffuse")
    plt.plot(theta_r_in_degrees_array, [x[0] + x[1] for x in fitted_data], label="fitted intensity")

    plt.title(title)
    plt.legend(bbox_to_anchor=(0.9, 0.9), bbox_transform=plt.gcf().transFigure)
    plt.ylabel("relative intensity")
    plt.xlabel("viewing angle (degrees)")
    string = "theta_i: " + str(theta_i_in_degrees) + "\nrho_L: " + str(rho_L) + "\nn: " + \
             str(n) + "\ngamma: " + str(gamma)
    plt.annotate(string, xy=(0.05, 0.8), xycoords='axes fraction')

plt.show()
