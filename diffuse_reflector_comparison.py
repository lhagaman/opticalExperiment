

import numpy as np
import matplotlib.pyplot as plt
from Point import Point
import TSTR_fit
from plotting import plot_with_TSTR_fit, make_points, make_data_by_run, plot_with_TSTR_fit_and_fitted_angles, plot_points


# solid angle of photodiode:
# 5.435 inches from sample
# 9 mm diameter
# in inches
distance_from_sample_to_photodiode = 5.435

# mm / (mm / inch)
photodiode_radius = (9 / 2.0) / 25.4

photodiode_solid_angle = np.pi * np.power(photodiode_radius, 2) / np.power(distance_from_sample_to_photodiode, 2)

photodiode_angular_width = photodiode_radius / distance_from_sample_to_photodiode

make_all_points = True
if make_all_points:
    # volts * amps/volt
    flux_i = 0.005830 * 100e-6

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

    data_coarse_30 = \
    make_data_by_run("10_18_diffuse_reflectors/diffuse reflector 120 grit in mineral oil 30 45 60.txt", -90, 90,
                     intensity_factor)[2]
    data_coarse_45 = \
    make_data_by_run("10_18_diffuse_reflectors/diffuse reflector 120 grit in mineral oil 30 45 60.txt", -90, 90,
                     intensity_factor)[1]
    data_coarse_60 = \
    make_data_by_run("10_18_diffuse_reflectors/diffuse reflector 120 grit in mineral oil 30 45 60.txt", -90, 90,
                     intensity_factor)[0]

    points_coarse_30 = make_points(data_coarse_30[0], 0, 30, 1.47399, 0.5, data_coarse_30[1], 405,
                                   photodiode_solid_angle, photodiode_angular_width, "30 degrees, 120 grit sample, older data")
    points_coarse_45 = make_points(data_coarse_45[0], 0, 45, 1.47399, 0.5, data_coarse_45[1], 405,
                                   photodiode_solid_angle, photodiode_angular_width, "45 degrees, 120 grit sample, older data")
    points_coarse_60 = make_points(data_coarse_60[0], 0, 60, 1.47399, 0.5, data_coarse_60[1], 405,
                                   photodiode_solid_angle, photodiode_angular_width, "60 degrees, 120 grit sample, older data")

    # volts * amps/volt
    flux_i = 0.00594 * 100e-6

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

    data_1_30 = \
    make_data_by_run("11_13_diffuse_reflectors/120 grit 30 45 60 75 11_13.txt", -90, 90,
                     intensity_factor)[3]
    data_1_45 = \
    make_data_by_run("11_13_diffuse_reflectors/120 grit 30 45 60 75 11_13.txt", -90, 90,
                     intensity_factor)[2]
    data_1_60 = \
    make_data_by_run("11_13_diffuse_reflectors/120 grit 30 45 60 75 11_13.txt", -90, 90,
                     intensity_factor)[1]
    data_1_75 = \
    make_data_by_run("11_13_diffuse_reflectors/120 grit 30 45 60 75 11_13.txt", -90, 90,
                     intensity_factor)[0]

    points_1_30 = make_points(data_1_30[0], 0, 30, 1.47399, 0.5, data_1_30[1], 405,
                                 photodiode_solid_angle, photodiode_angular_width, "30 degrees, 120 grit sample")
    points_1_45 = make_points(data_1_45[0], 0, 45, 1.47399, 0.5, data_1_45[1], 405,
                                 photodiode_solid_angle, photodiode_angular_width, "45 degrees, 120 grit sample")
    points_1_60 = make_points(data_1_60[0], 0, 60, 1.47399, 0.5, data_1_60[1], 405,
                                 photodiode_solid_angle, photodiode_angular_width, "60 degrees, 120 grit sample")
    points_1_75 = make_points(data_1_75[0], 0, 75, 1.47399, 0.5, data_1_75[1], 405,
                                 photodiode_solid_angle, photodiode_angular_width, "75 degrees, 120 grit sample")

    data_2_30 = \
        make_data_by_run("11_13_diffuse_reflectors/120 grit rotated 45 deg 30 45 60 75 11_13.txt", -90, 90,
                         intensity_factor)[3]
    data_2_45 = \
        make_data_by_run("11_13_diffuse_reflectors/120 grit rotated 45 deg 30 45 60 75 11_13.txt", -90, 90,
                         intensity_factor)[2]
    data_2_60 = \
        make_data_by_run("11_13_diffuse_reflectors/120 grit rotated 45 deg 30 45 60 75 11_13.txt", -90, 90,
                         intensity_factor)[1]
    data_2_75 = \
        make_data_by_run("11_13_diffuse_reflectors/120 grit rotated 45 deg 30 45 60 75 11_13.txt", -90, 90,
                         intensity_factor)[0]

    points_2_30 = make_points(data_2_30[0], 0, 30, 1.47399, 0.5, data_2_30[1], 405,
                              photodiode_solid_angle, photodiode_angular_width,
                              "30 degrees, 120 grit sample, rotated 45 degrees")
    points_2_45 = make_points(data_2_45[0], 0, 45, 1.47399, 0.5, data_2_45[1], 405,
                              photodiode_solid_angle, photodiode_angular_width,
                              "45 degrees, 120 grit sample, rotated 45 degrees")
    points_2_60 = make_points(data_2_60[0], 0, 60, 1.47399, 0.5, data_2_60[1], 405,
                              photodiode_solid_angle, photodiode_angular_width,
                              "60 degrees, 120 grit sample, rotated 45 degrees")
    points_2_75 = make_points(data_2_75[0], 0, 75, 1.47399, 0.5, data_2_75[1], 405,
                              photodiode_solid_angle, photodiode_angular_width,
                              "75 degrees, 120 grit sample, rotated 45 degrees")

all_points_coarse = points_coarse_30 + points_coarse_45 + points_coarse_60

all_points_1 = points_1_30 + points_1_45 + points_1_60 + points_1_75

all_points_2 = points_2_30 + points_2_45 + points_2_60 + points_2_75

all_points = all_points_coarse + all_points_1 + all_points_2

make_plots = False
if make_plots:
    plot_points(points_coarse_30 + points_1_30 + points_2_30, "120 Grit Diffuse Reflector 30 Degrees", log=False, show=False)
    plot_points(points_coarse_45 + points_1_45 + points_2_45, "120 Grit Diffuse Reflector 45 Degrees", log=False, show=False)
    plot_points(points_coarse_60 + points_1_60 + points_2_60, "120 Grit Diffuse Reflector 60 Degrees", log=False, show=False)
    plot_points(points_1_75 + points_2_75, "120 Grit Diffuse Reflector 75 Degrees", log=False, show=False)
    plt.show()



"""
points_coarse is the older data
points_1 is the data before rotating
points_2 is the data after rotating

"""

example_x_data = [p.theta_r_in_degrees for p in points_1_45]
example_y_data = [p.intensity for p in points_1_45]
plt.plot(example_x_data, example_y_data)
plt.show()


