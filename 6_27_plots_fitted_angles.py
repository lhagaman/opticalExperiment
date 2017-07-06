# currently not working at all, I think the function is too complicated for scipy to fit with these 4 variables

import numpy as np
import matplotlib.pyplot as plt
from Point import Point
import TSTR_fit


def make_points(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity_array,
                 wavelength, photodiode_solid_angle, run_name):
    points = []
    numpoints = len(theta_r_in_degrees_array)
    for i in range(numpoints):
        theta_r_in_degrees = theta_r_in_degrees_array[i]
        intensity = intensity_array[i]
        point = Point(theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity,
                 wavelength, photodiode_solid_angle, run_name)
        points.append(point)
    return points

# solid angle of photodiode:
# 5.435 inches from sample
# 9 mm diameter
# in inches
distance_from_sample_to_photodiode = 5.435

# mm / (mm / inch)
photodiode_radius = (9 / 2.0) / 25.4

photodiode_solid_angle = np.pi * np.power(photodiode_radius, 2) / np.power(distance_from_sample_to_photodiode, 2)

air_in_cylinder = True
if air_in_cylinder:

    # volts * amps/volt
    flux_i = 0.00186 * 100e-6

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

    dir = "cylindrical_cell_no_polarizer_6-28/"

    data_30 = np.loadtxt(dir + "30 deg.txt", skiprows=1)
    data_45 = np.loadtxt(dir + "45 deg.txt", skiprows=1)
    data_55 = np.loadtxt(dir + "55 deg.txt", skiprows=1)

    data_30_x = [np.round(point[0], 2) for point in data_30]
    data_30_y = [intensity_factor * point[1] for point in data_30]

    data_45_x = [np.round(point[0], 2) for point in data_45]
    data_45_y = [intensity_factor * point[1] for point in data_45]

    data_55_x = [np.round(point[0], 2) for point in data_55]
    data_55_y = [intensity_factor * point[1] for point in data_55]

    plot_data = False

    if plot_data:
        plt.plot(data_30_x, data_30_y, label="30 degrees")
        plt.plot(data_45_x, data_45_y, label="45 degrees")
        plt.plot(data_55_x, data_55_y, label="55 degrees")

        plt.title("Air, with tube, no polarizer")
        plt.legend()

        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.xlabel("viewing angle (degrees)")
        plt.show()

    points_30 = make_points(data_30_x, 0, 30, 1, 0.5, data_30_y, 405, photodiode_solid_angle,
                                     "30 degrees")
    points_45 = make_points(data_45_x, 0, 45, 1, 0.5, data_45_y, 405, photodiode_solid_angle,
                                     "45 degrees")
    points_55 = make_points(data_55_x, 0, 60, 1, 0.5, data_55_y, 405, photodiode_solid_angle,
                                     "55 degrees")

    points = points_30 + points_45 + points_55

    # cutting off regions which are dominated by reflections other than sample
    # off glass for low angles, and by background reflection for angles past parallel to sample surface
    cutoff_points_30 = []
    cutoff_data_30_x = []
    cutoff_data_30_y = []
    for point in points_30:
        if 10 < point.theta_r_in_degrees < 80:
            cutoff_points_30.append(point)
            cutoff_data_30_x.append(point.theta_r_in_degrees)
            cutoff_data_30_y.append(point.intensity)

    cutoff_points_45 = []
    cutoff_data_45_x = []
    cutoff_data_45_y = []
    for point in points_45:
        if 0 < point.theta_r_in_degrees < 80:
            cutoff_points_45.append(point)
            cutoff_data_45_x.append(point.theta_r_in_degrees)
            cutoff_data_45_y.append(point.intensity)

    cutoff_points_55 = []
    cutoff_data_55_x = []
    cutoff_data_55_y = []
    for point in points_55:
        if -20 < point.theta_r_in_degrees < 80:
            cutoff_points_55.append(point)
            cutoff_data_55_x.append(point.theta_r_in_degrees)
            cutoff_data_55_y.append(point.intensity)

    all_cutoff_points = cutoff_points_30 + cutoff_points_45 + cutoff_points_55

    adjusted_cutoff_points_30 = []
    for point in cutoff_points_30:
        adjusted_cutoff_points_30.append(Point(point.theta_r_in_degrees, 0, 35, 1, 0.5, point.intensity,
                     point.wavelength, point.photodiode_solid_angle, point.run_name))

    adjusted_cutoff_points_45 = []
    for point in cutoff_points_45:
        adjusted_cutoff_points_45.append(Point(point.theta_r_in_degrees, 0, 49, 1, 0.5, point.intensity,
                     point.wavelength, point.photodiode_solid_angle, point.run_name))

    adjusted_cutoff_points_55 = []
    for point in cutoff_points_55:
        adjusted_cutoff_points_55.append(Point(point.theta_r_in_degrees, 0, 57, 1, 0.5, point.intensity,
                     point.wavelength, point.photodiode_solid_angle, point.run_name))

    TSTR_parameters_30 = TSTR_fit.fit_parameters_and_angle(adjusted_cutoff_points_30)
    TSTR_parameters_45 = TSTR_fit.fit_parameters_and_angle(adjusted_cutoff_points_45)
    TSTR_parameters_55 = TSTR_fit.fit_parameters_and_angle(adjusted_cutoff_points_55)

    fit_30 = TSTR_fit.BRIDF_plotter(cutoff_data_30_x, 0, TSTR_parameters_30[0], 1, 0.5, TSTR_parameters_30[1:])
    fit_45 = TSTR_fit.BRIDF_plotter(cutoff_data_45_x, 0, TSTR_parameters_45[0], 1, 0.5, TSTR_parameters_45[1:])
    fit_55 = TSTR_fit.BRIDF_plotter(cutoff_data_55_x, 0, TSTR_parameters_55[0], 1, 0.5, TSTR_parameters_55[1:])

    plt.scatter(cutoff_data_30_x, cutoff_data_30_y, s=5)
    plt.scatter(cutoff_data_45_x, cutoff_data_45_y, s=5)
    plt.scatter(cutoff_data_55_x, cutoff_data_55_y, s=5)

    plt.semilogy(cutoff_data_30_x, fit_30)
    plt.semilogy(cutoff_data_45_x, fit_45)
    plt.semilogy(cutoff_data_55_x, fit_55)

    string = ""

    string += "theta_i: " + str(TSTR_parameters_30[0]) + "\n"
    string += "TSTR Parameters:\n"
    string += "rho_L: " + str(TSTR_parameters_30[1]) + ", n: " + str(TSTR_parameters_30[2]) + ", gamma: " + str(TSTR_parameters_30[3]) + "\n\n"

    string += "theta_i: " + str(TSTR_parameters_45[0]) + "\n"
    string += "TSTR Parameters:\n"
    string += "rho_L: " + str(TSTR_parameters_45[1]) + ", n: " + str(TSTR_parameters_45[2]) + ", gamma: " + str(TSTR_parameters_45[3]) + "\n\n"

    string += "theta_i: " + str(TSTR_parameters_55[0]) + "\n"
    string += "TSTR Parameters:\n"
    string += "rho_L: " + str(TSTR_parameters_55[1]) + ", n: " + str(TSTR_parameters_55[2]) + ", gamma: " + str(TSTR_parameters_55[3])

    plt.legend()
    plt.xlabel("viewing angle (degrees)")
    plt.ylabel("intensity (flux/str)/(input flux)")
    plt.annotate(string, xy=(0.05, 0.2), xycoords='axes fraction', size=8)
    plt.title("cylindrical cell in air, \nfitted incident angles, each run has own fit parameters")
    plt.show()

water = True
if water:

    # volts * amps/volt
    # estimated, changed for different angles of mirror in tube
    flux_i = 0.004 * 50e-6

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

    dir = "water_6-28/"

    data_30 = np.loadtxt(dir + "30 deg.txt", skiprows=1)
    data_45 = np.loadtxt(dir + "45 deg.txt", skiprows=1)
    data_55 = np.loadtxt(dir + "55 deg.txt", skiprows=1)

    data_30_x = [np.round(point[0], 2) for point in data_30]
    data_30_y = [intensity_factor * point[1] for point in data_30]

    data_45_x = [np.round(point[0], 2) for point in data_45]
    data_45_y = [intensity_factor * point[1] for point in data_45]

    data_55_x = [np.round(point[0], 2) for point in data_55]
    data_55_y = [intensity_factor * point[1] for point in data_55]

    plot_data = False

    if plot_data:
        plt.plot(data_30_x, data_30_y, label="30 degrees")
        plt.plot(data_45_x, data_45_y, label="45 degrees")
        plt.plot(data_55_x, data_55_y, label="55 degrees")

        plt.title("water, no polarizer")
        plt.legend()

        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.xlabel("viewing angle (degrees)")
        plt.show()

    points_30 = make_points(data_30_x, 0, 30, 1.33, 0.5, data_30_y, 405, photodiode_solid_angle,
                            "30 degrees")
    points_45 = make_points(data_45_x, 0, 45, 1.33, 0.5, data_45_y, 405, photodiode_solid_angle,
                            "45 degrees")
    points_55 = make_points(data_55_x, 0, 55, 1.33, 0.5, data_55_y, 405, photodiode_solid_angle,
                            "55 degrees")

    points = points_30 + points_45 + points_55

    # cutting off regions which are dominated by reflections other than sample
    # off glass for low angles, and by background reflection for angles past parallel to sample surface
    cutoff_points_30 = []
    cutoff_data_30_x = []
    cutoff_data_30_y = []
    for point in points_30:
        if 10 < point.theta_r_in_degrees < 80:
            cutoff_points_30.append(point)
            cutoff_data_30_x.append(point.theta_r_in_degrees)
            cutoff_data_30_y.append(point.intensity)

    cutoff_points_45 = []
    cutoff_data_45_x = []
    cutoff_data_45_y = []
    for point in points_45:
        if 0 < point.theta_r_in_degrees < 80:
            cutoff_points_45.append(point)
            cutoff_data_45_x.append(point.theta_r_in_degrees)
            cutoff_data_45_y.append(point.intensity)

    cutoff_points_55 = []
    cutoff_data_55_x = []
    cutoff_data_55_y = []
    for point in points_55:
        if -20 < point.theta_r_in_degrees < 80:
            cutoff_points_55.append(point)
            cutoff_data_55_x.append(point.theta_r_in_degrees)
            cutoff_data_55_y.append(point.intensity)

    all_cutoff_points = cutoff_points_30 + cutoff_points_45 + cutoff_points_55

    adjusted_cutoff_points_30 = []
    for point in cutoff_points_30:
        adjusted_cutoff_points_30.append(Point(point.theta_r_in_degrees, 0, 30, 1.33, 0.5, point.intensity,
                                               point.wavelength, point.photodiode_solid_angle, point.run_name))

    adjusted_cutoff_points_45 = []
    for point in cutoff_points_45:
        adjusted_cutoff_points_45.append(Point(point.theta_r_in_degrees, 0, 45, 1.33, 0.5, point.intensity,
                                               point.wavelength, point.photodiode_solid_angle, point.run_name))

    adjusted_cutoff_points_55 = []
    for point in cutoff_points_55:
        adjusted_cutoff_points_55.append(Point(point.theta_r_in_degrees, 0, 55, 1.33, 0.5, point.intensity,
                                               point.wavelength, point.photodiode_solid_angle, point.run_name))

    TSTR_parameters_30 = TSTR_fit.fit_parameters(adjusted_cutoff_points_30)
    TSTR_parameters_45 = TSTR_fit.fit_parameters(adjusted_cutoff_points_45)
    TSTR_parameters_55 = TSTR_fit.fit_parameters(adjusted_cutoff_points_55)

    fit_30 = TSTR_fit.BRIDF_plotter(cutoff_data_30_x, 0, 30, 1.33, 0.5, TSTR_parameters_30)
    fit_45 = TSTR_fit.BRIDF_plotter(cutoff_data_45_x, 0, 45, 1.33, 0.5, TSTR_parameters_45)
    fit_60 = TSTR_fit.BRIDF_plotter(cutoff_data_55_x, 0, 55, 1.33, 0.5, TSTR_parameters_55)

    plt.scatter(cutoff_data_30_x, cutoff_data_30_y, s=5)
    plt.scatter(cutoff_data_45_x, cutoff_data_45_y, s=5)
    plt.scatter(cutoff_data_55_x, cutoff_data_55_y, s=5)

    plt.semilogy(cutoff_data_30_x, fit_30)
    plt.semilogy(cutoff_data_45_x, fit_45)
    plt.semilogy(cutoff_data_55_x, fit_60)

    string = ""

    string += "theta_i: 30\n"
    string += "TSTR Parameters:\n"
    string += "rho_L: " + str(TSTR_parameters_30[0]) + ", n: " + str(TSTR_parameters_30[1]) + ", gamma: " + str(
        TSTR_parameters_30[2]) + "\n\n"

    string += "theta_i: 45\n"
    string += "TSTR Parameters:\n"
    string += "rho_L: " + str(TSTR_parameters_45[0]) + ", n: " + str(TSTR_parameters_45[1]) + ", gamma: " + str(
        TSTR_parameters_45[2]) + "\n\n"

    string += "theta_i: 55\n"
    string += "TSTR Parameters:\n"
    string += "rho_L: " + str(TSTR_parameters_55[0]) + ", n: " + str(TSTR_parameters_55[1]) + ", gamma: " + str(
        TSTR_parameters_55[2])

    plt.legend()
    plt.xlabel("viewing angle (degrees)")
    plt.ylabel("intensity (flux/str)/(input flux)")
    plt.annotate(string, xy=(0.05, 0.2), xycoords='axes fraction', size=8)
    plt.title(
        "cylindrical cell in water, \neach run has own fit parameters")
    plt.show()

mineral_oil = True
if mineral_oil:

    # volts * amps/volt
    # estimated, changed for different angles of mirror in tube
    flux_i = 0.0045 * 50e-6

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

    dir = "mineral_oil_6-28/"

    data_30 = np.loadtxt(dir + "30 deg.txt", skiprows=1)
    data_45 = np.loadtxt(dir + "45 deg.txt", skiprows=1)
    data_55 = np.loadtxt(dir + "55 deg.txt", skiprows=1)

    data_30_x = [np.round(point[0], 2) for point in data_30]
    data_30_y = [intensity_factor * point[1] for point in data_30]

    data_45_x = [np.round(point[0], 2) for point in data_45]
    data_45_y = [intensity_factor * point[1] for point in data_45]

    data_55_x = [np.round(point[0], 2) for point in data_55]
    data_55_y = [intensity_factor * point[1] for point in data_55]

    plot_data = False

    if plot_data:
        plt.plot(data_30_x, data_30_y, label="30 degrees")
        plt.plot(data_45_x, data_45_y, label="45 degrees")
        plt.plot(data_55_x, data_55_y, label="55 degrees")

        plt.title("mineral oil")
        plt.legend()

        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.xlabel("viewing angle (degrees)")
        plt.show()

    points_30 = make_points(data_30_x, 0, 30, 1.461, 0.5, data_30_y, 405, photodiode_solid_angle,
                            "30 degrees")
    points_45 = make_points(data_45_x, 0, 45, 1.461, 0.5, data_45_y, 405, photodiode_solid_angle,
                            "45 degrees")
    points_55 = make_points(data_55_x, 0, 65, 1.461, 0.5, data_55_y, 405, photodiode_solid_angle,
                            "55 degrees")

    points = points_30 + points_45 + points_55

    # cutting off regions which are dominated by reflections other than sample
    # off glass for low angles, and by background reflection for angles past parallel to sample surface
    cutoff_points_30 = []
    cutoff_data_30_x = []
    cutoff_data_30_y = []
    for point in points_30:
        if 10 < point.theta_r_in_degrees < 80:
            cutoff_points_30.append(point)
            cutoff_data_30_x.append(point.theta_r_in_degrees)
            cutoff_data_30_y.append(point.intensity)

    cutoff_points_45 = []
    cutoff_data_45_x = []
    cutoff_data_45_y = []
    for point in points_45:
        if 0 < point.theta_r_in_degrees < 80:
            cutoff_points_45.append(point)
            cutoff_data_45_x.append(point.theta_r_in_degrees)
            cutoff_data_45_y.append(point.intensity)

    cutoff_points_55 = []
    cutoff_data_55_x = []
    cutoff_data_55_y = []
    for point in points_55:
        if -20 < point.theta_r_in_degrees < 80:
            cutoff_points_55.append(point)
            cutoff_data_55_x.append(point.theta_r_in_degrees)
            cutoff_data_55_y.append(point.intensity)

    all_cutoff_points = cutoff_points_30 + cutoff_points_45 + cutoff_points_55

    adjusted_cutoff_points_30 = []
    for point in cutoff_points_30:
        adjusted_cutoff_points_30.append(Point(point.theta_r_in_degrees, 0, 30, 1.461, 0.5, point.intensity,
                                               point.wavelength, point.photodiode_solid_angle, point.run_name))

    adjusted_cutoff_points_45 = []
    for point in cutoff_points_45:
        adjusted_cutoff_points_45.append(Point(point.theta_r_in_degrees, 0, 45, 1.461, 0.5, point.intensity,
                                               point.wavelength, point.photodiode_solid_angle, point.run_name))

    adjusted_cutoff_points_55 = []
    for point in cutoff_points_55:
        adjusted_cutoff_points_55.append(Point(point.theta_r_in_degrees, 0, 69, 1.461, 0.5, point.intensity,
                                               point.wavelength, point.photodiode_solid_angle, point.run_name))

    TSTR_parameters_30 = TSTR_fit.fit_parameters(adjusted_cutoff_points_30)
    TSTR_parameters_45 = TSTR_fit.fit_parameters(adjusted_cutoff_points_45)
    TSTR_parameters_55 = TSTR_fit.fit_parameters(adjusted_cutoff_points_55)

    fit_30 = TSTR_fit.BRIDF_plotter(cutoff_data_30_x, 0, 35, 1.461, 0.5, TSTR_parameters_30)
    fit_45 = TSTR_fit.BRIDF_plotter(cutoff_data_45_x, 0, 45, 1.461, 0.5, TSTR_parameters_45)
    fit_60 = TSTR_fit.BRIDF_plotter(cutoff_data_55_x, 0, 69, 1.461, 0.5, TSTR_parameters_55)

    plt.scatter(cutoff_data_30_x, cutoff_data_30_y, s=5)
    plt.scatter(cutoff_data_45_x, cutoff_data_45_y, s=5)
    plt.scatter(cutoff_data_55_x, cutoff_data_55_y, s=5)

    plt.semilogy(cutoff_data_30_x, fit_30)
    plt.semilogy(cutoff_data_45_x, fit_45)
    plt.semilogy(cutoff_data_55_x, fit_60)

    string = ""

    string += "theta_i: 30\n"
    string += "TSTR Parameters:\n"
    string += "rho_L: " + str(TSTR_parameters_30[0]) + ", n: " + str(TSTR_parameters_30[1]) + ", gamma: " + str(
        TSTR_parameters_30[2]) + "\n\n"

    string += "theta_i: 45\n"
    string += "TSTR Parameters:\n"
    string += "rho_L: " + str(TSTR_parameters_45[0]) + ", n: " + str(TSTR_parameters_45[1]) + ", gamma: " + str(
        TSTR_parameters_45[2]) + "\n\n"

    string += "theta_i: 69\n"
    string += "TSTR Parameters:\n"
    string += "rho_L: " + str(TSTR_parameters_55[0]) + ", n: " + str(TSTR_parameters_55[1]) + ", gamma: " + str(
        TSTR_parameters_55[2])

    plt.legend()
    plt.xlabel("viewing angle (degrees)")
    plt.ylabel("intensity (flux/str)/(input flux)")
    plt.annotate(string, xy=(0.05, 0.2), xycoords='axes fraction', size=8)
    plt.title(
        "cylindrical cell in mineral oil, \neach run has own fit parameters")
    plt.show()

