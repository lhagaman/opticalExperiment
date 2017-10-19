

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
    flux_i_vertical = 0.002093 * 100e-6
    flux_i_horizontal = 0.003907 * 100e-6

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor_vertical = sensitivity / (photodiode_solid_angle * flux_i_vertical)
    intensity_factor_horizontal = sensitivity / (photodiode_solid_angle * flux_i_horizontal)

    data_45_vertical = make_data_by_run("10_12_polarization_data/45 degree vertical and horizontal mineral oil.txt", -90, 90, intensity_factor_vertical)[1]
    data_45_horizontal = make_data_by_run("10_12_polarization_data/45 degree vertical and horizontal mineral oil.txt", -90, 90, intensity_factor_horizontal)[0]
    data_60_horizontal = make_data_by_run("10_12_polarization_data/60 degree horizontal and vertical mineral oil.txt", -90, 90, intensity_factor_horizontal)[1]
    data_60_vertical = make_data_by_run("10_12_polarization_data/60 degree horizontal and vertical mineral oil.txt", -90, 90, intensity_factor_vertical)[0]
    data_75_vertical = make_data_by_run("10_12_polarization_data/75 degree vertical and horizontal mineral oil.txt", -90, 90, intensity_factor_vertical)[1]
    data_75_horizontal = make_data_by_run("10_12_polarization_data/75 degree vertical and horizontal mineral oil.txt", -90, 90, intensity_factor_horizontal)[0]

    points_45_vertical = make_points(data_45_vertical[0], 0, 45, 1.47399, 0.5, data_45_vertical[1], 405,
                                     photodiode_solid_angle, photodiode_angular_width, "45 degrees vertical polarization")
    points_45_horizontal = make_points(data_45_horizontal[0], 0, 45, 1.47399, 0.5, data_45_horizontal[1], 405,
                                     photodiode_solid_angle, photodiode_angular_width, "45 degrees horizontal polarization")
    points_60_vertical = make_points(data_60_vertical[0], 0, 60, 1.47399, 0.5, data_60_vertical[1], 405,
                                     photodiode_solid_angle, photodiode_angular_width,
                                     "60 degrees vertical polarization")
    points_60_horizontal = make_points(data_60_horizontal[0], 0, 60, 1.47399, 0.5, data_60_horizontal[1], 405,
                                       photodiode_solid_angle, photodiode_angular_width,
                                       "60 degrees horizontal polarization")
    points_75_vertical = make_points(data_75_vertical[0], 0, 75, 1.47399, 0.5, data_75_vertical[1], 405,
                                     photodiode_solid_angle, photodiode_angular_width,
                                     "75 degrees vertical polarization")
    points_75_horizontal = make_points(data_75_horizontal[0], 0, 75, 1.47399, 0.5, data_75_horizontal[1], 405,
                                       photodiode_solid_angle, photodiode_angular_width,
                                       "75 degrees horizontal polarization")
all_points = points_45_vertical + points_45_horizontal + points_60_vertical + points_60_horizontal + points_75_vertical + points_75_horizontal

plot_points(points_45_horizontal + points_45_vertical, "45 degrees", log=False, show=False)
plot_points(points_60_horizontal + points_60_vertical, "60 degrees", log=False, show=False)
plot_points(points_75_horizontal + points_75_vertical, "75 degrees", log=False, show=False)
plt.show()

