import TSTR_fit
import gaussian_fit
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


plot_TSTR_unpolarized = False
plot_TSTR_polarizations = False
plot_gaussian = False
plot_TSTR_unpolarized_and_gaussian = True

if plot_TSTR_unpolarized:

    plt.figure()

    title = "45 degree IBA in air, Torrance-Sparrow-Trowbridge-Reitz Fit, unpolarized"

    polarization = 0

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

    parameters = TSTR_fit.fit_parameters(points)

    fitted_total_data = TSTR_fit.BRIDF_plotter(
        theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, parameters)
    fitted_specular_data = TSTR_fit.BRIDF_specular_plotter(
        theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, parameters)
    fitted_diffuse_data = TSTR_fit.BRIDF_diffuse_plotter(
        theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, parameters)

    plt.plot(theta_r_in_degrees_array, intensity_array, label="experimental")
    plt.plot(theta_r_in_degrees_array, fitted_specular_data, label="fitted specular")
    plt.plot(theta_r_in_degrees_array, fitted_diffuse_data, label="fitted diffuse")
    plt.plot(theta_r_in_degrees_array, fitted_total_data, label="fitted intensity")

    plt.title(title)
    plt.legend(bbox_to_anchor=(0.9, 0.9), bbox_transform=plt.gcf().transFigure)
    plt.ylabel("relative intensity")
    plt.xlabel("viewing angle (degrees)")
    rho_L = parameters[0]
    n = parameters[1]
    gamma = parameters[2]
    string = "theta_i: " + str(theta_i_in_degrees) + "\nrho_L: " + str(rho_L) + "\nn: " + \
             str(n) + "\ngamma: " + str(gamma)
    plt.annotate(string, xy=(0.05, 0.8), xycoords='axes fraction')

    plt.show()

if plot_TSTR_polarizations:

    for polarization in [0, 1, 2]:

        plt.figure(polarization)

        if (polarization == 0):
            p = "unpolarized"
        elif (polarization == 1):
            p = "s-polarized"
        elif (polarization == 2):
            p = "p-polarized"
        title = "45 degree IBA in air, TSTR Fit, " + p

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

        parameters = TSTR_fit.fit_parameters(points)

        rho_L = parameters[0]
        n = parameters[1]
        gamma = parameters[2]

        fitted_data = TSTR_fit.BRIDF_plotter(
            theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, parameters)

        plt.plot(theta_r_in_degrees_array, intensity_array, label="experimental")
        plt.plot(theta_r_in_degrees_array, fitted_data, label="fitted intensity")

        plt.title(title)
        plt.legend(bbox_to_anchor=(0.9, 0.9), bbox_transform=plt.gcf().transFigure)
        plt.ylabel("relative intensity")
        plt.xlabel("viewing angle (degrees)")
        string = "theta_i: " + str(theta_i_in_degrees) + "\nrho_L: " + str(rho_L) + "\nn: " + \
                 str(n) + "\ngamma: " + str(gamma)
        plt.annotate(string, xy=(0.05, 0.8), xycoords='axes fraction')

if plot_gaussian:

    plt.figure()

    # always arbitrary for gaussian fit
    # (actual polarization will likely affect R1 and R2 though)
    polarization = 0

    title = "45 degree IBA in air, Gaussian Fit"

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

    parameters = gaussian_fit.fit_parameters(points)

    fitted_total_data = gaussian_fit.BRIDF_plotter(theta_r_in_degrees_array, theta_i_in_degrees, parameters)
    fitted_specular_data = gaussian_fit.BRIDF_specular_plotter(theta_r_in_degrees_array,
                                                               theta_i_in_degrees, parameters)
    fitted_diffuse_data = gaussian_fit.BRIDF_diffuse_plotter(theta_r_in_degrees_array, theta_i_in_degrees, parameters)

    plt.plot(theta_r_in_degrees_array, intensity_array, label="experimental")
    plt.plot(theta_r_in_degrees_array, fitted_specular_data, label="fitted specular")
    plt.plot(theta_r_in_degrees_array, fitted_diffuse_data, label="fitted diffuse")
    plt.plot(theta_r_in_degrees_array, fitted_total_data, label="fitted intensity")

    plt.title(title)
    plt.legend(bbox_to_anchor=(0.9, 0.9), bbox_transform=plt.gcf().transFigure)
    plt.ylabel("relative intensity")
    plt.xlabel("viewing angle (degrees)")

    sigma = parameters[0]
    R_1 = parameters[1]
    R_2 = parameters[2]

    string = "theta_i: " + str(theta_i_in_degrees) + "\nsigma: " + str(sigma) + "\nR_1: " + \
             str(R_1) + "\nR_2: " + str(R_2)
    plt.annotate(string, xy=(0.05, 0.8), xycoords='axes fraction')

if plot_TSTR_unpolarized_and_gaussian:

    title = "45 degree IBA in air, TSTR unpolarized and Gaussian Fit"

    polarization = 0

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

    gaussian_parameters = gaussian_fit.fit_parameters(points)
    sigma = gaussian_parameters[0]
    R_1 = gaussian_parameters[1]
    R_2 = gaussian_parameters[2]

    TSTR_parameters = TSTR_fit.fit_parameters(points)
    rho_L = TSTR_parameters[0]
    n = TSTR_parameters[1]
    gamma = TSTR_parameters[2]

    gaussian_fitted_data = gaussian_fit.BRIDF_plotter(theta_r_in_degrees_array, theta_i_in_degrees, gaussian_parameters)

    TSTR_fitted_data = TSTR_fit.BRIDF_plotter(theta_r_in_degrees_array, phi_r_in_degrees,
                                        theta_i_in_degrees, n_0, polarization, TSTR_parameters)

    plt.plot(theta_r_in_degrees_array, intensity_array, label="experimental")
    plt.plot(theta_r_in_degrees_array, gaussian_fitted_data, label="Gaussian model")
    plt.plot(theta_r_in_degrees_array, TSTR_fitted_data,
             label="TSTR model")

    plt.title(title)
    plt.legend(bbox_to_anchor=(0.9, 0.9), bbox_transform=plt.gcf().transFigure)
    plt.ylabel("relative intensity")
    plt.xlabel("viewing angle (degrees)")
    string = "theta_i: " + str(theta_i_in_degrees) + "\n\nsigma: " + str(sigma) + "\nR_1: " + \
             str(R_1) + "\nR_2: " + str(R_2) + "\n\nrho_L: " + str(rho_L) + "\nn: " + \
                 str(n) + "\ngamma: " + str(gamma)
    plt.annotate(string, xy=(0.05, 0.6), xycoords='axes fraction')

plt.show()
