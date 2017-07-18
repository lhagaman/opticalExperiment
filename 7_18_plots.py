
# air through tube, 2.05 mV, 100uA/V
# water through tube, 2.12 mV, 100 uA/V
# mineral oil through tube: 2.38 mV, 100uA/V


import numpy as np
import matplotlib.pyplot as plt
from Point import Point
import TSTR_fit
from plotting import plot_with_TSTR_fit, make_points, make_data_by_run, plot_with_TSTR_fit_and_fitted_angles

# in inches
distance_from_sample_to_photodiode = 5.435

# these are guessed, need to measure later
photodiode_height = 0.450
photodiode_width = 0.075

photodiode_area = photodiode_height * photodiode_width

photodiode_solid_angle = photodiode_area / np.power(distance_from_sample_to_photodiode, 2)

make_all_points = True
if make_all_points:
    make_air_points = True
    if make_air_points:
        # volts * amps/volt
        flux_i = 0.00205 * 100e-6

        # amps per volt during measurements
        sensitivity = 100 * 1e-9

        # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
        # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
        intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

        data_30 = make_data_by_run("7_18_plots/30 deg in air with slit.txt", 10, 85, intensity_factor)[0]
        data_45 = make_data_by_run("7_18_plots/45 deg in air with slit.txt", 20, 80, intensity_factor)[0]
        data_60 = make_data_by_run("7_18_plots/60 in air with slit.txt", 0, 80, intensity_factor)[0]
        data_75 = make_data_by_run("7_18_plots/75 in air with slit.txt", 0, 80, intensity_factor)[0]

        points_30 = make_points(data_30[0], 0, 30, 1, 0.5, data_30[1], 405, photodiode_solid_angle,
                                "30 degrees air with slit")
        points_45 = make_points(data_45[0], 0, 45, 1, 0.5, data_45[1], 405, photodiode_solid_angle,
                                "45 degrees air with slit")
        points_60 = make_points(data_60[0], 0, 60, 1, 0.5, data_60[1], 405, photodiode_solid_angle,
                                "60 degrees air with slit")
        points_75 = make_points(data_75[0], 0, 75, 1, 0.5, data_75[1], 405, photodiode_solid_angle,
                                "75 degrees air slit")

        points = points_30 + points_45 + points_60 + points_75
        points_no_75 = points_30 + points_45 + points_60

fit_plot_air = True
if fit_plot_air:
    plot_with_TSTR_fit_and_fitted_angles(points, "TSTR Fit in Air With Slit")
