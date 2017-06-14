from plotting import plot_with_semi_empirical_and_gaussian_fits, plot_semi_empirical_components, \
    plot_large_gas_layer_gaussian
import numpy as np
import matplotlib.pyplot as plt
from Point import Point


read_data = True
plot_all_original_data = False

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


# I assume perpendicular beam to windows, only perfectly true for cylindrical vessel
def frac_reflected_at_window(n_0, n_1, n_2):
    frac_reflected_first_boundary = np.power((n_0 - n_1) / (n_0 + n_1), 2)
    frac_reflected_second_boundary = np.power((n_1 - n_2) / (n_1 + n_2), 2)
    # I only include the two shortest paths for a photon to get reflected, could make more precise
    return frac_reflected_first_boundary + (1 - frac_reflected_first_boundary) * frac_reflected_second_boundary * \
                                           (1 - frac_reflected_first_boundary)


frac_reflected_air = frac_reflected_at_window(air_index_of_refraction, CaF2_index_of_refraction,
                                              air_index_of_refraction) * \
                     frac_reflected_at_window(air_index_of_refraction, CaF2_index_of_refraction,
                                              air_index_of_refraction)

frac_reflected_water = frac_reflected_at_window(air_index_of_refraction, CaF2_index_of_refraction,
                                                water_index_of_refraction) * \
                       frac_reflected_at_window(water_index_of_refraction, CaF2_index_of_refraction,
                                                air_index_of_refraction)

# going through two windows
frac_transmitted_air = 0.94 * 0.94

# frac_absorbed is purely a property of CaF2, should be independent of medium
# frac_transmitted = frac_not_reflected * frac_not_absorbed
frac_absorbed = -frac_transmitted_air / (1 - frac_reflected_air) + 1

frac_transmitted_water = (1 - frac_absorbed) * (1 - frac_reflected_water)

flux_i_air = frac_transmitted_air * frac_transmitted_air * flux_i_without_windows

flux_i_water = frac_transmitted_water * frac_transmitted_water * flux_i_without_windows

# amps per volt during measurements
sensitivity = 100 * 1e-9

# product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
# intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
air_intensity_factor = sensitivity / (photodiode_solid_angle * flux_i_air)

water_intensity_factor = intensity_factor = sensitivity / (photodiode_solid_angle * flux_i_water)

# CaF2 viewport Lesker, look at transmittivity, take that into account

if read_data:

    dir = "405 nm Blue Laser diode reflectance measurements/"

    data_30_degree_1 = np.loadtxt(dir + "30 deg with blue laser diode run 2 after angle change 6-6.txt", skiprows=1)

    data_30_degree_2 = np.loadtxt(dir + "30 deg with blue laser diode run 2 after angle change.txt", skiprows=1)

    data_30_degree_3 = np.loadtxt(dir + "30 deg with blue laser diode.txt", skiprows=1)

    data_45_degree_1 = np.loadtxt(dir + "45 deg with blue laser diode 6-5 run 2 after changing angles.txt", skiprows=1)

    data_45_degree_2 = np.loadtxt(dir + "45 deg with blue laser diode 6-5.txt", skiprows=1)

    data_45_degree_3 = np.loadtxt(dir + "45 deg with blue laser diode.txt", skiprows=1)

    data_45_degree_water = np.loadtxt(dir + "45 deg with blue laser diode in water  6-6.txt", skiprows=1)

    data_30_degree_1_x = [np.round(point[0], 2) for point in data_30_degree_1]
    data_30_degree_1_y = [intensity_factor * point[1] for point in data_30_degree_1]

    data_30_degree_2_x = [np.round(point[0], 2) for point in data_30_degree_2]
    data_30_degree_2_y = [intensity_factor * point[1] for point in data_30_degree_2]

    data_30_degree_3_x = [np.round(point[0], 2) for point in data_30_degree_3]
    data_30_degree_3_y = [intensity_factor * point[1] for point in data_30_degree_3]

    data_45_degree_1_x = [np.round(point[0], 2) for point in data_45_degree_1]
    data_45_degree_1_y = [intensity_factor * point[1] for point in data_45_degree_1]

    data_45_degree_2_x = [np.round(point[0], 2) for point in data_45_degree_2]
    data_45_degree_2_y = [intensity_factor * point[1] for point in data_45_degree_2]

    data_45_degree_3_x = [np.round(point[0], 2) for point in data_45_degree_3]
    data_45_degree_3_y = [intensity_factor * point[1] for point in data_45_degree_3]

    data_45_degree_water_x = [np.round(point[0], 2) for point in data_45_degree_water]
    data_45_degree_water_y = [water_intensity_factor * point[1] for point in data_45_degree_water]


if plot_all_original_data:

    f, axarr = plt.subplots(2, sharex=False)

    axarr[0].plot(data_30_degree_1_x, data_30_degree_1_y, label="run 2 after angle change 6-6")
    axarr[0].plot(data_30_degree_2_x, data_30_degree_2_y, label="run 2 after angle change")
    axarr[0].plot(data_30_degree_3_x, data_30_degree_3_y, label="before angle change")

    axarr[0].set_title("30 deg with blue laser diode")
    axarr[0].legend()

    axarr[1].plot(data_45_degree_1_x, data_45_degree_1_y, label="6-5 run 2 after changing angles")
    axarr[1].plot(data_45_degree_2_x, data_45_degree_2_y, label="6-5")
    axarr[1].plot(data_45_degree_3_x, data_45_degree_3_y, label="run 1")
    axarr[1].plot(data_45_degree_water_x, data_45_degree_water_y, label="in water")

    axarr[1].set_title("45 deg with blue laser diode")
    axarr[1].legend()

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

points_30_degree_1 = make_points(data_30_degree_1_x, 0, 30, 1, 0.5, data_30_degree_1_y, 405, photodiode_solid_angle,
                                 "30 deg with blue laser diode run 2 after angle change 6-6")
points_30_degree_2 = make_points(data_30_degree_2_x, 0, 30, 1, 0.5, data_30_degree_2_y, 405, photodiode_solid_angle,
                                 "30 deg with blue laser diode run 2 after angle change")
points_30_degree_3 = make_points(data_30_degree_3_x, 0, 30, 1, 0.5, data_30_degree_3_y, 405, photodiode_solid_angle,
                                 "30 deg with blue laser diode.txt")
points_45_degree_1 = make_points(data_45_degree_1_x, 0, 45, 1, 0.5, data_45_degree_1_y, 405, photodiode_solid_angle,
                                 "45 deg with blue laser diode 6-5 run 2 after changing angles")
points_45_degree_2 = make_points(data_45_degree_2_x, 0, 45, 1, 0.5, data_45_degree_2_y, 405, photodiode_solid_angle,
                                 "45 deg with blue laser diode 6-5")
points_45_degree_3 = make_points(data_45_degree_3_x, 0, 45, 1, 0.5, data_45_degree_3_y, 405, photodiode_solid_angle,
                                 "45 deg with blue laser diode")
points_45_degree_water = make_points(data_45_degree_water_x, 0, 45, 1.33, 0.5, data_45_degree_water_y, 405,
                                     photodiode_solid_angle, "45 deg with blue laser diode in water  6-6")

all_points = points_30_degree_1 + points_30_degree_2 + points_30_degree_3 + points_45_degree_1 + \
             points_45_degree_2 + points_45_degree_3

# here I manually adjust some input angles to see how it affects the fit
# because I suspect some have large errors
adjusted_points_30_degree_1 = make_points(data_30_degree_1_x, 0, 32, 1, 0.5, data_30_degree_1_y, 405, photodiode_solid_angle,
                                 "30 deg with blue laser diode run 2 after angle change 6-6")
adjusted_points_30_degree_2 = make_points(data_30_degree_2_x, 0, 32, 1, 0.5, data_30_degree_2_y, 405, photodiode_solid_angle,
                                 "30 deg with blue laser diode run 2 after angle change")
adjusted_points_30_degree_3 = make_points(data_30_degree_3_x, 0, 41, 1, 0.5, data_30_degree_3_y, 405, photodiode_solid_angle,
                                 "30 deg with blue laser diode.txt")
adjusted_points_45_degree_1 = make_points(data_45_degree_1_x, 0, 46, 1, 0.5, data_45_degree_1_y, 405, photodiode_solid_angle,
                                 "45 deg with blue laser diode 6-5 run 2 after changing angles")
adjusted_points_45_degree_2 = make_points(data_45_degree_2_x, 0, 50, 1, 0.5, data_45_degree_2_y, 405, photodiode_solid_angle,
                                 "45 deg with blue laser diode 6-5")
adjusted_points_45_degree_3 = make_points(data_45_degree_3_x, 0, 52, 1, 0.5, data_45_degree_3_y, 405, photodiode_solid_angle,
                                 "45 deg with blue laser diode")

all_points_adjusted = adjusted_points_30_degree_1 + adjusted_points_30_degree_2 + adjusted_points_30_degree_3 + \
                      adjusted_points_45_degree_1 + adjusted_points_45_degree_2 + adjusted_points_45_degree_3

# here I cut off small and large angles where I suspect interference from the edges of the window
cutoff_45_degree_points_adjusted = []
for point in all_points_adjusted:
    if 30 <= point.theta_r_in_degrees <= 70 and 43 <= point.theta_i_in_degrees <= 55:
        cutoff_45_degree_points_adjusted.append(point)

cutoff_adjusted_45_degree_1 = []
for point in adjusted_points_45_degree_1:
    if 30 <= point.theta_r_in_degrees <= 70:
        cutoff_adjusted_45_degree_1.append(point)

cutoff_adjusted_45_degree_2 = []
for point in adjusted_points_45_degree_1:
    if 30 <= point.theta_r_in_degrees <= 70:
        cutoff_adjusted_45_degree_2.append(point)

cutoff_adjusted_45_degree_3 = []
for point in adjusted_points_45_degree_1:
    if 30 <= point.theta_r_in_degrees <= 70:
        cutoff_adjusted_45_degree_3.append(point)

cutoff_45_degree_water = []
for point in points_45_degree_water:
    if 30 <= point.theta_r_in_degrees <= 70:
        cutoff_45_degree_water.append(point)

# these fit poorly
# plot_with_semi_empirical_and_gaussian_fits(points_30_degree_1,
#                                           "30 deg with blue laser diode run 2 after angle change 6-6.txt")
# plot_with_semi_empirical_and_gaussian_fits(points_45_degree_2, "45 deg with blue laser diode 6-5")
# plot_with_semi_empirical_and_gaussian_fits(points_45_degree_3, "45 deg with blue laser diode")

# these fit well
# plot_with_semi_empirical_and_gaussian_fits(points_45_degree_1)
# plot_with_semi_empirical_and_gaussian_fits(adjusted_points_45_degree_2,
#                                           "45 deg with blue laser diode 6-5, adjusted theta_i")
# plot_with_semi_empirical_and_gaussian_fits(adjusted_points_45_degree_3,
#                                           "45 deg with blue laser diode, adjusted theta_i")
# plot_with_semi_empirical_and_gaussian_fits(adjusted_points_45_degree_4, "45 deg, adjusted theta_i")

# plot_with_semi_empirical_and_gaussian_fits(all_points, "All Points")
# plot_with_semi_empirical_and_gaussian_fits(all_points_adjusted, "All Points With Manually Adjusted theta_i Values")

# plot_with_semi_empirical_and_gaussian_fits(cutoff_45_degree_points_adjusted,
#                                           "45 degree points, cutoff edge angles, manually adjusted theta_i values")

# plot_with_semi_empirical_and_gaussian_fits(cutoff_adjusted_45_degree_1, "45 deg with blue laser diode 6-5 run 2 "
#                                                                        "after changing angles, \ncutoff edge angles, "
#                                                                        "manually adjusted theta_i value")
# plot_with_semi_empirical_and_gaussian_fits(cutoff_adjusted_45_degree_2)
# plot_with_semi_empirical_and_gaussian_fits(cutoff_adjusted_45_degree_3)

# plot_with_semi_empirical_and_gaussian_fits(cutoff_45_degree_water)


# points = cutoff_adjusted_45_degree_1 + cutoff_45_degree_water
# plot_with_semi_empirical_and_gaussian_fits(points)

# plot_semi_empirical_components(cutoff_adjusted_45_degree_1)
# plot_semi_empirical_components(cutoff_45_degree_water)

plot_large_gas_layer_gaussian(cutoff_adjusted_45_degree_1)

