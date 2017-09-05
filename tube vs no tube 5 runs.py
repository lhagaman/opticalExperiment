# sample 1 VP through tube, 0.0005240 V
# sample 1 HP through tube, 0.002138 V
# sample 2 VP through tube, 0.0005660 V
# sample 2 HP through tube, 0.002195 V
# sample 2 HP no tube: 0.002819 V

import numpy as np
import matplotlib.pyplot as plt
import TSTR_fit
from plotting import plot_with_TSTR_fit, make_points, plot_TSTR_fit_one_set_of_parameters, make_points_std, make_data_by_run, make_data_all_runs, plot_with_TSTR_fit_and_fitted_angles, plot_points, plot_with_semi_empirical_fit, plot_points_error_bars, plot_TSTR_fit_error_bars_no_show

# in inches
distance_from_sample_to_photodiode = 5.435

photodiode_height = 0.450
photodiode_width = 0.075

photodiode_angular_width_radians = 0.075 / 5.435
photodiode_angular_width = photodiode_angular_width_radians * 180. / np.pi

photodiode_area = photodiode_height * photodiode_width

photodiode_solid_angle = photodiode_area / np.power(distance_from_sample_to_photodiode, 2)

make_all_points = True
if make_all_points:

    theta_i = 45

    # these correct for the polarizer and polarization of the laser
    flux_i_air_horizontal_sample_2 = 0.002195 * 100e-6
    flux_i_air_horizontal_sample_2_no_tube = 0.002819 * 100e-6

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i

    intensity_factor_air_horizontal_sample_2 = sensitivity / (photodiode_solid_angle * flux_i_air_horizontal_sample_2)
    intensity_factor_air_horizontal_sample_2_no_tube = sensitivity / (photodiode_solid_angle * flux_i_air_horizontal_sample_2_no_tube)

    data_45_horizontal_air_sample_2 = make_data_all_runs("five runs/Air/45 deg HP new sample.txt", 0, 85, intensity_factor_air_horizontal_sample_2)
    points_45_horizontal_air_sample_2 = make_points_std(data_45_horizontal_air_sample_2[0], 0, theta_i, 1, 1, data_45_horizontal_air_sample_2[1], 405, photodiode_solid_angle, photodiode_angular_width,
                                                    "horizontal sample 2")

    data_45_horizontal_air_sample_2_no_tube = make_data_all_runs("five runs/Air/45 deg HP new sample no tube.txt", 0, 85, intensity_factor_air_horizontal_sample_2_no_tube)
    points_45_horizontal_air_sample_2_no_tube = make_points_std(data_45_horizontal_air_sample_2_no_tube[0], 0, theta_i, 1, 1, data_45_horizontal_air_sample_2_no_tube[1], 405, photodiode_solid_angle,
                                                        photodiode_angular_width, "horizontal sample 2 no tube")

    all_points = points_45_horizontal_air_sample_2 + points_45_horizontal_air_sample_2_no_tube

plot_1 = True
if plot_1:
    plot_points_error_bars(all_points, "With vs. Without Cylinder, Air, 45 Degrees, Horizontal Polarization")
