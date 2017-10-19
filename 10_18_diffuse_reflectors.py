

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

    data_fine_30 = make_data_by_run("10_18_diffuse_reflectors/diffuse reflector 1500 grit in mineral oil 30, 45, 60, 75.txt", -90, 90, intensity_factor)[3]
    data_fine_45 = make_data_by_run("10_18_diffuse_reflectors/diffuse reflector 1500 grit in mineral oil 30, 45, 60, 75.txt", -90, 90, intensity_factor)[2]
    data_fine_60 = make_data_by_run("10_18_diffuse_reflectors/diffuse reflector 1500 grit in mineral oil 30, 45, 60, 75.txt", -90, 90, intensity_factor)[1]
    data_fine_75 = make_data_by_run("10_18_diffuse_reflectors/diffuse reflector 1500 grit in mineral oil 75 full range.txt", -90, 90, intensity_factor)[0]

    data_coarse_30 = make_data_by_run("10_18_diffuse_reflectors/diffuse reflector 120 grit in mineral oil 30 45 60.txt", -90, 90, intensity_factor)[2]
    data_coarse_45 = make_data_by_run("10_18_diffuse_reflectors/diffuse reflector 120 grit in mineral oil 30 45 60.txt", -90, 90, intensity_factor)[1]
    data_coarse_60 = make_data_by_run("10_18_diffuse_reflectors/diffuse reflector 120 grit in mineral oil 30 45 60.txt", -90, 90, intensity_factor)[0]

    points_fine_30 = make_points(data_fine_30[0], 0, 30, 1.47399, 0.5, data_fine_30[1], 405,
                                 photodiode_solid_angle, photodiode_angular_width, "30 degrees, 1500 grit sample")
    points_fine_45 = make_points(data_fine_45[0], 0, 45, 1.47399, 0.5, data_fine_45[1], 405,
                                 photodiode_solid_angle, photodiode_angular_width, "45 degrees, 1500 grit sample")
    points_fine_60 = make_points(data_fine_60[0], 0, 60, 1.47399, 0.5, data_fine_60[1], 405,
                                 photodiode_solid_angle, photodiode_angular_width, "60 degrees, 1500 grit sample")
    points_fine_75 = make_points(data_fine_75[0], 0, 75, 1.47399, 0.5, data_fine_75[1], 405,
                                 photodiode_solid_angle, photodiode_angular_width, "75 degrees, 1500 grit sample")

    points_coarse_30 = make_points(data_coarse_30[0], 0, 30, 1.47399, 0.5, data_coarse_30[1], 405,
                                 photodiode_solid_angle, photodiode_angular_width, "30 degrees, 120 grit sample")
    points_coarse_45 = make_points(data_coarse_45[0], 0, 45, 1.47399, 0.5, data_coarse_45[1], 405,
                                 photodiode_solid_angle, photodiode_angular_width, "45 degrees, 120 grit sample")
    points_coarse_60 = make_points(data_coarse_60[0], 0, 60, 1.47399, 0.5, data_coarse_60[1], 405,
                                 photodiode_solid_angle, photodiode_angular_width, "60 degrees, 120 grit sample")

all_points_fine = points_fine_30 + points_fine_45 + points_fine_60 + points_fine_75

all_points_coarse = points_coarse_30 + points_coarse_45 + points_coarse_60

plot_points(all_points_fine, "1500 Grit Diffuse Reflector", log=False, show=False)
plot_points(all_points_coarse, "120 Grit Diffuse Reflector", log=False, show=False)
plt.show()
