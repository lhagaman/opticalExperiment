# tested with filter over 405 nm laser

from plotting import plot_with_semi_empirical_and_gaussian_fits, plot_semi_empirical_components, \
    plot_large_gas_layer_gaussian, plot_partial_gas_layer_gaussian, \
    plot_with_semi_empirical_TSTR_gaussian_and_partial_gas_layer_fits, plot_with_reff_polynomial_fit
import numpy as np
import matplotlib.pyplot as plt
from Point import Point

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

# in amps measured from photodiode (0.533 volts at 2 microamps / volt)
# measured through chopper wheel, not through windows
flux_i_without_windows = 0.533 * 2e-6

# ~94% transmission for CaF2 windows according to
# https://www.lesker.com/newweb/flanges/viewports_cf_exotic.cfm?pgid=0
# for non-air medium, this should be more carefully calculated

# according to http://www.sydor.com/wp-content/uploads/Corning-Calcium-Fluoride-CaF2.pdf
CaF2_index_of_refraction = 1.44
water_index_of_refraction = 1.33
air_index_of_refraction = 1

# amps per volt during measurements
sensitivity = 10e-9

# laser power in amps from photodiode
laser_power_through_windows = 380e-12

# product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
# intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
intensity_factor = sensitivity / (photodiode_solid_angle * laser_power_through_windows)

dir = ""

data = np.loadtxt(dir + "preamp sens 10 na-V 45 deg and background runs 150 hz.txt", skiprows=1)

rounded_data = []
for point in data:
    rounded_data.append([np.round(point[0], 2), intensity_factor * point[1]])

data_grouped_by_run_number = [[], [], [], []]
i = -1
for point in rounded_data:
    if point[0] == 16:
        i += 1
    data_grouped_by_run_number[i].append(point)

background_1 = data_grouped_by_run_number[0]
background_2 = data_grouped_by_run_number[1]
data_1 = data_grouped_by_run_number[2]
data_2 = data_grouped_by_run_number[3]

background_1_x = [point[0] for point in background_1]
background_1_y = [point[1] for point in background_1]

background_2_x = [point[0] for point in background_2]
background_2_y = [point[1] for point in background_2]

data_1_x = [point[0] for point in data_1]
data_1_y = [point[1] for point in data_1]

data_2_x = [point[0] for point in data_2]
data_2_y = [point[1] for point in data_2]

background_average_x = background_1_x
background_average_y = [(background_1_y[i] + background_2_y[i]) / 2.0 for i in range(len(background_1_y))]

data_average_x = data_1_x
data_average_y = [(data_1_y[i] + data_2_y[i]) / 2.0 for i in range(len(data_1_y))]

subtracted_data_x = data_average_x
subtracted_data_y = [data_average_y[i] - background_average_y[i] for i in range(len(data_1_y))]


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

subtracted_points = make_points(subtracted_data_x, 0, 45, 1, 0.5, subtracted_data_y, 405, photodiode_solid_angle,
                                 "45 deg 405 nm with filter")

"""
plot_with_semi_empirical_TSTR_gaussian_and_partial_gas_layer_fits(subtracted_points)

plt.plot(background_average_x, background_average_y, label="background")
plt.plot(data_average_x, data_average_y, label="data")
plt.legend()
plt.show()
"""

plot_with_reff_polynomial_fit(subtracted_points)
