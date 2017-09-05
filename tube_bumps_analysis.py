# sample 1 VP through tube, 0.0005240 V
# sample 1 HP through tube, 0.002138 V
# sample 2 VP through tube, 0.0005660 V
# sample 2 HP through tube, 0.002195

import numpy as np
import matplotlib.pyplot as plt
import TSTR_fit
from plotting import plot_with_TSTR_fit, make_points, plot_TSTR_fit_one_set_of_parameters, make_points_std, make_data_by_run, make_data_by_run_angles_count_down, make_data_all_runs, plot_with_TSTR_fit_and_fitted_angles, plot_points, plot_with_semi_empirical_fit, plot_points_error_bars, plot_TSTR_fit_error_bars_no_show

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

    # these are from a different day, we didn't measure it for this test, we average them to compare
    flux_i_air_vertical_sample_1 = 0.0005240 * 100e-6
    flux_i_air_horizontal_sample_1 = 0.002138 * 100e-6

    flux_i = (flux_i_air_horizontal_sample_1 + flux_i_air_vertical_sample_1) / 2.

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

    data_75 = make_data_by_run("tube bumps analysis data/relative to normal angles.txt", -90, 90, intensity_factor)[0]
    data_45 = make_data_by_run("tube bumps analysis data/relative to normal angles.txt", -90, 90, intensity_factor)[1]
    data_60 = make_data_by_run("tube bumps analysis data/relative to normal angles.txt", -90, 90, intensity_factor)[2]
    data_30 = make_data_by_run("tube bumps analysis data/relative to normal angles.txt", -90, 90, intensity_factor)[3]

    points_30 = make_points(data_30[0], 0, 30, 1, 0.5, data_30[1], 405, photodiode_solid_angle, photodiode_angular_width, "30 degrees")
    points_45 = make_points(data_45[0], 0, 45, 1, 0.5, data_45[1], 405, photodiode_solid_angle, photodiode_angular_width, "45 degrees")
    points_60 = make_points(data_60[0], 0, 60, 1, 0.5, data_60[1], 405, photodiode_solid_angle, photodiode_angular_width, "60 degrees")
    points_75 = make_points(data_75[0], 0, 75, 1, 0.5, data_75[1], 405, photodiode_solid_angle, photodiode_angular_width, "75 degrees")

    data_75_rot = make_data_by_run_angles_count_down("tube bumps analysis data/rotation stage angles.txt", 0, 180, intensity_factor)[0]
    data_45_rot = make_data_by_run_angles_count_down("tube bumps analysis data/rotation stage angles.txt", 0, 180, intensity_factor)[1]
    data_60_rot = make_data_by_run_angles_count_down("tube bumps analysis data/rotation stage angles.txt", 0, 180, intensity_factor)[2]
    data_30_rot = make_data_by_run_angles_count_down("tube bumps analysis data/rotation stage angles.txt", 0, 180, intensity_factor)[3]

    points_30_rot = make_points(data_30_rot[0], 0, 30, 1, 0.5, data_30_rot[1], 405, photodiode_solid_angle, photodiode_angular_width, "30 degrees, relative to rotation stage")
    points_45_rot = make_points(data_45_rot[0], 0, 45, 1, 0.5, data_45_rot[1], 405, photodiode_solid_angle, photodiode_angular_width, "45 degrees, relative to rotation stage")
    points_60_rot = make_points(data_60_rot[0], 0, 60, 1, 0.5, data_60_rot[1], 405, photodiode_solid_angle, photodiode_angular_width, "60 degrees, relative to rotation stage")
    points_75_rot = make_points(data_75_rot[0], 0, 75, 1, 0.5, data_75_rot[1], 405, photodiode_solid_angle, photodiode_angular_width, "75 degrees, relative to rotation stage")

    all_points = points_30 + points_45 + points_60 + points_75

    all_points_rot = points_30_rot + points_45_rot + points_60_rot + points_75_rot

plot_1 = False
if plot_1:
    #plot_points(points_30_rot, "30 degrees")
    plot_points(all_points_rot, "all points rotation stage angle")
    plot_points(all_points, "all points")
    plt.show()

fit = True
if fit:
    fit_30 = TSTR_fit.fit_parameters(points_30)
    fit_45 = TSTR_fit.fit_parameters(points_45)
    fit_60 = TSTR_fit.fit_parameters(points_60)
    fit_75 = TSTR_fit.fit_parameters(points_75)

plot_1b = False
if plot_1b:
    plot_TSTR_fit_one_set_of_parameters(points_30, fit_30, "30 degrees")
    plot_TSTR_fit_one_set_of_parameters(points_45, fit_45, "45 degrees")
    plot_TSTR_fit_one_set_of_parameters(points_60, fit_60, "60 degrees")
    plot_TSTR_fit_one_set_of_parameters(points_75, fit_75, "75 degrees")


plot_2 = True
if plot_2:
    points_30_x = [point.theta_r_in_degrees for point in points_30]
    points_30_y = [point.intensity for point in points_30]

    points_30_x_rot = [point.theta_r_in_degrees for point in points_30_rot]
    points_30_y_rot = [point.intensity for point in points_30]

    fit_30_x = points_30_x
    fit_30_y = TSTR_fit.BRIDF_plotter(fit_30_x, 0, 30, 1, 0.5, fit_30)
    fit_30_x_rot = [180. - (30 + x) for x in fit_30_x]
    fit_30_y_rot = fit_30_y

    differences_30 = [points_30_y_rot[i] - fit_30_y_rot[i] for i in range(len(points_30_y_rot))]
    ratios_30 = [differences_30[i] / fit_30_y_rot[i] for i in range(len(differences_30))]


    points_45_x = [point.theta_r_in_degrees for point in points_45]
    points_45_y = [point.intensity for point in points_45]

    points_45_x_rot = [point.theta_r_in_degrees for point in points_45_rot]
    points_45_y_rot = [point.intensity for point in points_45]

    fit_45_x = points_45_x
    fit_45_y = TSTR_fit.BRIDF_plotter(fit_45_x, 0, 45, 1, 0.5, fit_45)
    fit_45_x_rot = [180. - (45 + x) for x in fit_45_x]
    fit_45_y_rot = fit_45_y

    differences_45 = [points_45_y_rot[i] - fit_45_y_rot[i] for i in range(len(points_45_y_rot))]
    ratios_45 = [differences_45[i] / fit_45_y_rot[i] for i in range(len(differences_45))]


    points_60_x = [point.theta_r_in_degrees for point in points_60]
    points_60_y = [point.intensity for point in points_60]

    points_60_x_rot = [point.theta_r_in_degrees for point in points_60_rot]
    points_60_y_rot = [point.intensity for point in points_60]

    fit_60_x = points_60_x
    fit_60_y = TSTR_fit.BRIDF_plotter(fit_60_x, 0, 60, 1, 0.5, fit_60)
    fit_60_x_rot = [180. - (60 + x) for x in fit_60_x]
    fit_60_y_rot = fit_60_y

    differences_60 = [points_60_y_rot[i] - fit_60_y_rot[i] for i in range(len(points_60_y_rot))]
    ratios_60 = [differences_60[i] / fit_60_y_rot[i] for i in range(len(differences_60))]


    points_75_x = [point.theta_r_in_degrees for point in points_75]
    points_75_y = [point.intensity for point in points_75]

    points_75_x_rot = [point.theta_r_in_degrees for point in points_75_rot]
    points_75_y_rot = [point.intensity for point in points_75]

    fit_75_x = points_75_x
    fit_75_y = TSTR_fit.BRIDF_plotter(fit_75_x, 0, 75, 1, 0.5, fit_75)
    fit_75_x_rot = [180. - (75 + x) for x in fit_75_x]
    fit_75_y_rot = fit_75_y

    differences_75 = [points_75_y_rot[i] - fit_75_y_rot[i] for i in range(len(points_75_y_rot))]
    ratios_75 = [differences_75[i] / fit_75_y_rot[i] for i in range(len(differences_75))]


    plt.plot(points_30_x_rot, ratios_30, label="30 degrees", c="r")
    plt.plot(points_45_x_rot, ratios_45, label="45 degrees", c="g")
    plt.plot(points_60_x_rot, ratios_60, label="60 degrees", c="b")
    plt.plot(points_75_x_rot, ratios_75, label="75 degrees", c="m")

    plt.axvline(x=120, c="r")
    plt.axvline(x=90, c="g")
    plt.axvline(x=60, c="b")
    plt.axvline(x=30, c="m")


    plt.legend()
    plt.title("fractional difference from each fit\n(specular angles marked with vertical lines)")
    plt.xlabel("reflected angle measured relative to the rotation stage\n(0 means photodiode is opposite laser, 180 means back of photodiode is blocking  laser)")
    plt.ylabel("fractional difference between data and fit \n((data - fit)/fit)")
    plt.show()



