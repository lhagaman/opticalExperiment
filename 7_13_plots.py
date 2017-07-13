
# air through tube, 2.05 mV, 100uA/V
# water through tube, 2.12 mV, 100 uA/V
# mineral oil through tube: 2.38 mV, 100uA/V


import numpy as np
import matplotlib.pyplot as plt
from Point import Point
import TSTR_fit
from plotting import plot_with_TSTR_fit

# solid angle of photodiode:
# 5.435 inches from sample
# 9 mm diameter
# in inches
distance_from_sample_to_photodiode = 5.435

# mm / (mm / inch)
photodiode_radius = (9 / 2.0) / 25.4

photodiode_solid_angle = np.pi * np.power(photodiode_radius, 2) / np.power(distance_from_sample_to_photodiode, 2)


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

        dir = "7_13_data/"

        data = np.loadtxt(dir + "75, 60, 45, and 30 in air.txt", skiprows=1)
        data_by_run = []
        for i in range(len(data)):
            if i == 0 or data[i][0] < data[i - 1][0]:
                data_by_run.append([data[i]])
            else:
                data_by_run[-1].append(data[i])

        data_75 = data_by_run[3]
        data_60 = data_by_run[2]
        data_45 = data_by_run[1]
        data_30 = data_by_run[0]

        data_30_x = [np.round(point[0], 2) for point in data_30]
        data_30_y = [intensity_factor * point[1] for point in data_30]

        data_45_x = [np.round(point[0], 2) for point in data_45]
        data_45_y = [intensity_factor * point[1] for point in data_45]

        data_60_x = [np.round(point[0], 2) for point in data_60]
        data_60_y = [intensity_factor * point[1] for point in data_60]

        data_75_x = [np.round(point[0], 2) for point in data_75]
        data_75_y = [intensity_factor * point[1] for point in data_75]

        points_30 = make_points(data_30_x, 0, 30, 1, 0.5, data_30_y, 405, photodiode_solid_angle,
                                "30 degrees")
        points_45 = make_points(data_45_x, 0, 45, 1, 0.5, data_45_y, 405, photodiode_solid_angle,
                                "45 degrees")
        points_60 = make_points(data_60_x, 0, 60, 1, 0.5, data_60_y, 405, photodiode_solid_angle,
                                "60 degrees")
        points_75 = make_points(data_75_x, 0, 75, 1, 0.5, data_75_y, 405, photodiode_solid_angle,
                                "75 degrees")

        points = points_30 + points_45 + points_60 + points_75

        # cutting off regions which are dominated by reflections other than sample
        # off glass for low angles, and by background reflection for angles past parallel to sample surface
        cutoff_points_30 = []
        cutoff_data_30_x = []
        cutoff_data_30_y = []
        for point in points_30:
            if -10 < point.theta_r_in_degrees < 85:
                cutoff_points_30.append(point)
                cutoff_data_30_x.append(point.theta_r_in_degrees)
                cutoff_data_30_y.append(point.intensity)

        cutoff_points_45 = []
        cutoff_data_45_x = []
        cutoff_data_45_y = []
        for point in points_45:
            if -30 < point.theta_r_in_degrees < 85:
                cutoff_points_45.append(point)
                cutoff_data_45_x.append(point.theta_r_in_degrees)
                cutoff_data_45_y.append(point.intensity)

        cutoff_points_60 = []
        cutoff_data_60_x = []
        cutoff_data_60_y = []
        for point in points_60:
            if -30 < point.theta_r_in_degrees < 85:
                cutoff_points_60.append(point)
                cutoff_data_60_x.append(point.theta_r_in_degrees)
                cutoff_data_60_y.append(point.intensity)

        cutoff_points_75 = []
        cutoff_data_75_x = []
        cutoff_data_75_y = []
        for point in points_75:
            if -50 < point.theta_r_in_degrees < 85:
                cutoff_points_75.append(point)
                cutoff_data_75_x.append(point.theta_r_in_degrees)
                cutoff_data_75_y.append(point.intensity)

        all_cutoff_points_air = cutoff_points_30 + cutoff_points_45 + cutoff_points_60 + cutoff_points_75

    make_water_points = True
    if make_water_points:

        # volts * amps/volt
        flux_i = 0.00212 * 100e-6

        # amps per volt during measurements
        sensitivity = 100 * 1e-9

        # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
        # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
        intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

        dir = "7_13_data/"

        data = np.loadtxt(dir + "75, 60, 45, and 30 in water.txt", skiprows=1)
        data_by_run = []
        for i in range(len(data)):
            if i == 0 or data[i][0] < data[i - 1][0]:
                data_by_run.append([data[i]])
            else:
                data_by_run[-1].append(data[i])

        data_75 = data_by_run[3]
        data_60 = data_by_run[2]
        data_45 = data_by_run[1]
        data_30 = data_by_run[0]

        data_30_x = [np.round(point[0], 2) for point in data_30]
        data_30_y = [intensity_factor * point[1] for point in data_30]

        data_45_x = [np.round(point[0], 2) for point in data_45]
        data_45_y = [intensity_factor * point[1] for point in data_45]

        data_60_x = [np.round(point[0], 2) for point in data_60]
        data_60_y = [intensity_factor * point[1] for point in data_60]

        data_75_x = [np.round(point[0], 2) for point in data_75]
        data_75_y = [intensity_factor * point[1] for point in data_75]

        points_30 = make_points(data_30_x, 0, 30, 1.33, 0.5, data_30_y, 405, photodiode_solid_angle,
                                "30 degrees")
        points_45 = make_points(data_45_x, 0, 45, 1.33, 0.5, data_45_y, 405, photodiode_solid_angle,
                                "45 degrees")
        points_60 = make_points(data_60_x, 0, 60, 1.33, 0.5, data_60_y, 405, photodiode_solid_angle,
                                "60 degrees")
        points_75 = make_points(data_75_x, 0, 75, 1.33, 0.5, data_75_y, 405, photodiode_solid_angle,
                                "75 degrees")

        points = points_30 + points_45 + points_60 + points_75

        # cutting off regions which are dominated by reflections other than sample
        # off glass for low angles, and by background reflection for angles past parallel to sample surface
        cutoff_points_30 = []
        cutoff_data_30_x = []
        cutoff_data_30_y = []
        for point in points_30:
            if -10 < point.theta_r_in_degrees < 85:
                cutoff_points_30.append(point)
                cutoff_data_30_x.append(point.theta_r_in_degrees)
                cutoff_data_30_y.append(point.intensity)

        cutoff_points_45 = []
        cutoff_data_45_x = []
        cutoff_data_45_y = []
        for point in points_45:
            if -30 < point.theta_r_in_degrees < 85:
                cutoff_points_45.append(point)
                cutoff_data_45_x.append(point.theta_r_in_degrees)
                cutoff_data_45_y.append(point.intensity)

        cutoff_points_60 = []
        cutoff_data_60_x = []
        cutoff_data_60_y = []
        for point in points_60:
            if -30 < point.theta_r_in_degrees < 85:
                cutoff_points_60.append(point)
                cutoff_data_60_x.append(point.theta_r_in_degrees)
                cutoff_data_60_y.append(point.intensity)

        cutoff_points_75 = []
        cutoff_data_75_x = []
        cutoff_data_75_y = []
        for point in points_75:
            if -50 < point.theta_r_in_degrees < 85:
                cutoff_points_75.append(point)
                cutoff_data_75_x.append(point.theta_r_in_degrees)
                cutoff_data_75_y.append(point.intensity)

        all_cutoff_points_water = cutoff_points_30 + cutoff_points_45 + cutoff_points_60 + cutoff_points_75

    make_mineral_oil_points = True
    if make_mineral_oil_points:

        # volts * amps/volt
        flux_i = 0.00238 * 100e-6

        # amps per volt during measurements
        sensitivity = 100 * 1e-9

        # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
        # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
        intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

        dir = "7_13_data/"

        data = np.loadtxt(dir + "75, 60, 45, and 30 in mineral oil.txt", skiprows=1)
        data_by_run = []
        for i in range(len(data)):
            if i == 0 or data[i][0] < data[i - 1][0]:
                data_by_run.append([data[i]])
            else:
                data_by_run[-1].append(data[i])

        data_75 = data_by_run[3]
        data_60 = data_by_run[2]
        data_45 = data_by_run[1]
        data_30 = data_by_run[0]

        data_30_x = [np.round(point[0], 2) for point in data_30]
        data_30_y = [intensity_factor * point[1] for point in data_30]

        data_45_x = [np.round(point[0], 2) for point in data_45]
        data_45_y = [intensity_factor * point[1] for point in data_45]

        data_60_x = [np.round(point[0], 2) for point in data_60]
        data_60_y = [intensity_factor * point[1] for point in data_60]

        data_75_x = [np.round(point[0], 2) for point in data_75]
        data_75_y = [intensity_factor * point[1] for point in data_75]

        points_30 = make_points(data_30_x, 0, 30, 1.461, 0.5, data_30_y, 405, photodiode_solid_angle,
                                "30 degrees")
        points_45 = make_points(data_45_x, 0, 45, 1.461, 0.5, data_45_y, 405, photodiode_solid_angle,
                                "45 degrees")
        points_60 = make_points(data_60_x, 0, 60, 1.461, 0.5, data_60_y, 405, photodiode_solid_angle,
                                "60 degrees")
        points_75 = make_points(data_75_x, 0, 75, 1.461, 0.5, data_75_y, 405, photodiode_solid_angle,
                                "75 degrees")

        points = points_30 + points_45 + points_60 + points_75

        # cutting off regions which are dominated by reflections other than sample
        # off glass for low angles, and by background reflection for angles past parallel to sample surface
        cutoff_points_30 = []
        cutoff_data_30_x = []
        cutoff_data_30_y = []
        for point in points_30:
            if -10 < point.theta_r_in_degrees < 85:
                cutoff_points_30.append(point)
                cutoff_data_30_x.append(point.theta_r_in_degrees)
                cutoff_data_30_y.append(point.intensity)

        cutoff_points_45 = []
        cutoff_data_45_x = []
        cutoff_data_45_y = []
        for point in points_45:
            if -30 < point.theta_r_in_degrees < 85:
                cutoff_points_45.append(point)
                cutoff_data_45_x.append(point.theta_r_in_degrees)
                cutoff_data_45_y.append(point.intensity)

        cutoff_points_60 = []
        cutoff_data_60_x = []
        cutoff_data_60_y = []
        for point in points_60:
            if -30 < point.theta_r_in_degrees < 85:
                cutoff_points_60.append(point)
                cutoff_data_60_x.append(point.theta_r_in_degrees)
                cutoff_data_60_y.append(point.intensity)

        cutoff_points_75 = []
        cutoff_data_75_x = []
        cutoff_data_75_y = []
        for point in points_75:
            if -50 < point.theta_r_in_degrees < 85:
                cutoff_points_75.append(point)
                cutoff_data_75_x.append(point.theta_r_in_degrees)
                cutoff_data_75_y.append(point.intensity)

        all_cutoff_points_mineral_oil = cutoff_points_30 + cutoff_points_45 + cutoff_points_60 + cutoff_points_75

    all_cutoff_points = all_cutoff_points_air + all_cutoff_points_water + all_cutoff_points_mineral_oil

fit_plot_air = False
if fit_plot_air:
    plot_with_TSTR_fit(all_cutoff_points_air, "TSTR Fit in Air")

fit_plot_water = False
if fit_plot_water:
    plot_with_TSTR_fit(all_cutoff_points_water, "TSTR Fit in Water")

fit_plot_mineral_oil = False
if fit_plot_mineral_oil:
    plot_with_TSTR_fit(all_cutoff_points_mineral_oil, "TSTR Fit in Mineral Oil")

fit_plot_all_points = True
if fit_plot_all_points:
    plot_with_TSTR_fit(all_cutoff_points, "TSTR Fit in Air, Water, and Mineral Oil")
