# sample 1 VP through tube, 0.0005240 V
# sample 1 HP through tube, 0.002138 V
# sample 2 VP through tube, 0.0005660 V
# sample 2 HP through tube, 0.002195

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

    theta_i_1 = 47
    theta_i_2 = 50

    # these are all through the tube, they correct for the polarizer and polarization of the laser

    flux_i_air_vertical_sample_1 = 0.0005240 * 100e-6
    flux_i_air_horizontal_sample_1 = 0.002138 * 100e-6
    flux_i_air_vertical_sample_2 = 0.0005660 * 100e-6
    flux_i_air_horizontal_sample_2 = 0.002195 * 100e-6

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor_air_vertical_sample_1 = sensitivity / (photodiode_solid_angle * flux_i_air_vertical_sample_1)
    intensity_factor_air_horizontal_sample_1 = sensitivity / (photodiode_solid_angle * flux_i_air_horizontal_sample_1)

    intensity_factor_air_vertical_sample_2 = sensitivity / (photodiode_solid_angle * flux_i_air_vertical_sample_2)
    intensity_factor_air_horizontal_sample_2 = sensitivity / (photodiode_solid_angle * flux_i_air_horizontal_sample_2)

    data_45_horizontal_air_sample_1 = make_data_all_runs("five runs/Air/45 deg HP.txt", -90, 90, intensity_factor_air_horizontal_sample_1)
    points_45_horizontal_air_sample_1 = make_points_std(data_45_horizontal_air_sample_1[0], 0, theta_i_1, 1, 1, data_45_horizontal_air_sample_1[1], 405, photodiode_solid_angle, photodiode_angular_width,
                                                    "horizontal sample 1")

    data_45_horizontal_air_sample_2 = make_data_all_runs("five runs/Air/45 deg HP new sample.txt", -90, 90, intensity_factor_air_horizontal_sample_2)
    points_45_horizontal_air_sample_2 = make_points_std(data_45_horizontal_air_sample_2[0], 0, theta_i_2, 1, 1, data_45_horizontal_air_sample_2[1], 405, photodiode_solid_angle, photodiode_angular_width,
                                                    "horizontal sample 2")

    data_45_vertical_air_sample_1 = make_data_all_runs("five runs/Air/45 deg VP.txt", -90, 90, intensity_factor_air_vertical_sample_1)
    points_45_vertical_air_sample_1 = make_points_std(data_45_vertical_air_sample_1[0], 0, theta_i_1, 1, 0, data_45_vertical_air_sample_1[1], 405, photodiode_solid_angle, photodiode_angular_width,
                                                    "vertical sample 1")

    data_45_vertical_air_sample_2 = make_data_all_runs("five runs/Air/45 deg VP new sample.txt", -90, 90, intensity_factor_air_vertical_sample_2)
    points_45_vertical_air_sample_2 = make_points_std(data_45_vertical_air_sample_2[0], 0, theta_i_2, 1, 0, data_45_vertical_air_sample_2[1], 405, photodiode_solid_angle, photodiode_angular_width,
                                                    "vertical sample 2")

    all_points = points_45_vertical_air_sample_1 + points_45_horizontal_air_sample_1 + points_45_horizontal_air_sample_2 + points_45_vertical_air_sample_2

fit_1 = TSTR_fit.fit_std(points_45_horizontal_air_sample_1 + points_45_vertical_air_sample_1)
fit_2 = TSTR_fit.fit_std(points_45_horizontal_air_sample_2 + points_45_vertical_air_sample_2)


plot_1 = True
if plot_1:
    plot_TSTR_fit_error_bars_no_show(all_points, "Fits Of Different Polarizations For Two Different Samples, 45 Degrees, Air")
    plt.show()

# Ryan's plot
plot_2 = False
if plot_2:
    parameters = TSTR_fit.fit_std(points_45_vertical_air_sample_2)
    plot_TSTR_fit_one_set_of_parameters(points_45_vertical_air_sample_2 + points_45_horizontal_air_sample_2, parameters, "Sample 2")
    plt.show()



