from plotting import plot_with_semi_empirical_and_gaussian_fits, plot_semi_empirical_components, \
    plot_large_gas_layer_gaussian, plot_partial_gas_layer_gaussian, \
    plot_with_semi_empirical_TSTR_gaussian_and_partial_gas_layer_fits
import numpy as np
import matplotlib.pyplot as plt
from Point import Point
import TSTR_fit
import silva_fit

# solid angle of photodiode:
# 5.435 inches from sample
# 9 mm diameter
# total flux of laser: 0.533 volts * 2 microamps/volt
# measurements are multiplied by 100 nanoamps per volt

# in inches
distance_from_sample_to_photodiode = 5.435

# mm / (mm / inch)
photodiode_radius = (9 / 2.0) / 25.4

photodiode_solid_angle = np.pi * np.power(photodiode_radius, 2) / np.power(distance_from_sample_to_photodiode, 2)

# in amps measured from photodiode (0.152 millivots at 1 milliamp / volt)
# measured through chopper wheel, not through windows
flux_i = 0.152e-3 * 1e-3

# amps per volt during measurements
sensitivity = 100 * 1e-9

# product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
# intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

dir = "air_6-27/"

data_30 = np.loadtxt(dir + "30 deg in air.txt", skiprows=1)
data_45 = np.loadtxt(dir + "45 deg in air.txt", skiprows=1)
data_60 = np.loadtxt(dir + "60 deg in air.txt", skiprows=1)

data_30_x = [np.round(point[0], 2) for point in data_30]
data_30_y = [intensity_factor * point[1] for point in data_30]

data_45_x = [np.round(point[0], 2) for point in data_45]
data_45_y = [intensity_factor * point[1] for point in data_45]

data_60_x = [np.round(point[0], 2) for point in data_60]
data_60_y = [intensity_factor * point[1] for point in data_60]

plot_data = False

if plot_data:
    plt.plot(data_30_x, data_30_y, label="30 degrees")
    plt.plot(data_45_x, data_45_y, label="45 degrees")
    plt.plot(data_60_x, data_60_y, label="60 degrees")

    plt.title("Air, no tube, unknown linear polaization")
    plt.legend()

    plt.ylabel("intensity (flux/str)/(input flux)")
    plt.xlabel("viewing angle (degrees)")
    plt.show()


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

# assume unpolarized since polarization is unknown
points_30 = make_points(data_30_x, 0, 30, 1, 0.5, data_30_y, 405, photodiode_solid_angle,
                                 "30 degrees")
points_45 = make_points(data_30_x, 0, 45, 1, 0.5, data_30_y, 405, photodiode_solid_angle,
                                 "30 degrees")
points_60 = make_points(data_30_x, 0, 60, 1, 0.5, data_30_y, 405, photodiode_solid_angle,
                                 "30 degrees")

points = points_30 + points_45 + points_60


silva_parameters = silva_fit.fit_parameters(points)
"""
plt.scatter(data_30_x, data_30_y)
diffuse = []
lobe = []
spike = []
total = []
for point in points_30:
    trio = silva_fit.BRIDF_trio(point.theta_r_in_degrees * np.pi / 180, point.phi_r_in_degrees * np.pi / 180, point.theta_i_in_degrees * np.pi / 180, point.n_0, point.polarization, photodiode_solid_angle, silva_parameters)
    diffuse.append(trio[0])
    lobe.append(trio[1])
    spike.append(trio[2])
    total.append(trio[0] + trio[1] + trio[2])
plt.plot(data_30_x, diffuse, label="diffuse")
plt.plot(data_30_x, lobe, label="lobe")
plt.plot(data_30_x, spike, label="spike")
plt.plot(data_30_x, total, label="total")
plt.legend()
plt.show()
"""

diffuse = []
lobe = []
spike = []
total = []
point = points[0]
angles = [i * 10 for i in range(9)]

# from paper
# parameters has the form [rho_L, K, n, gamma]
parameters = [0.74, 1.45, 1.7, 0.049]
for theta_i_in_degrees in angles:
    print("theta_i: ", theta_i_in_degrees)
    diffuse.append(silva_fit.reflectance_diffuse(theta_i_in_degrees * np.pi / 180, point.n_0, point.polarization,
                                photodiode_solid_angle, parameters))


plt.plot(angles, diffuse, label="diffuse")
"""
plt.plot(angles, lobe, label="specular lobe")
plt.plot(angles, spike, label="specular spike")
plt.plot(angles, total, label="total")
"""
plt.xlabel("angle of incidence (degrees)")
plt.ylabel("reflectance")
plt.legend()
plt.show()


