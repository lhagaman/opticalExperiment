# air through tube, 2.05 mV, 100uA/V
# water through tube, 2.12 mV, 100 uA/V
# mineral oil through tube: 2.38 mV, 100uA/V


import numpy as np
import matplotlib.pyplot as plt
from Point import Point
import TSTR_fit
from plotting import plot_with_TSTR_fit, make_points, make_data_by_run, plot_with_TSTR_fit_and_fitted_angles, plot_points, plot_with_semi_empirical_fit

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

    # these are all through the tube, they correct for the polarizer and polarization of the laser

    flux_i_air_vertical = 0.0004890 * 100e-6
    flux_i_air_horizontal = 0.002120 * 100e-6

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor_air_vertical = sensitivity / (photodiode_solid_angle * flux_i_air_vertical)
    intensity_factor_air_horizontal = sensitivity / (photodiode_solid_angle * flux_i_air_horizontal)

    data_45_horizontal_air = make_data_by_run("HP and VP for 45 deg in air.txt", 0, 80, intensity_factor_air_horizontal)[1]

    points_45_horizontal_air = make_points(data_45_horizontal_air[0], 0, 47, 1, 1, data_45_horizontal_air[1], 405, photodiode_solid_angle, photodiode_angular_width,
                                             "45 degrees horizontal water")

    data_45_vertical_air = make_data_by_run("HP and VP for 45 deg in air.txt", 0, 80, intensity_factor_air_vertical)[0]
    points_45_vertical_air = make_points(data_45_vertical_air[0], 0, 47, 1, 0, data_45_vertical_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "45 degrees vertical air")

plot = True
if plot:
    plot_with_TSTR_fit(points_45_vertical_air + points_45_horizontal_air, "Different Polarizations in Air With Slit")
