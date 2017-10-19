

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
    flux_i_no_tube = 0.005850 * 100e-6
    flux_i_with_tube = 0.005940 * 100e-6
    flux_i_older = 0.002603 * 100e-6 # from 9/12/2017

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor_no_tube = sensitivity / (photodiode_solid_angle * flux_i_no_tube)
    intensity_factor_with_tube = sensitivity / (photodiode_solid_angle * flux_i_with_tube)
    intensity_factor_older = sensitivity / (photodiode_solid_angle * flux_i_older)


    data_old_walls = make_data_by_run("background in mineral oil no sample old walls new walls lens tubes.txt", -90, 90, intensity_factor_no_tube)[2]
    data_new_walls = make_data_by_run("background in mineral oil no sample old walls new walls lens tubes.txt", -90, 90, intensity_factor_no_tube)[1]
    data_with_tube = make_data_by_run("background in mineral oil no sample old walls new walls lens tubes.txt", -90, 90, intensity_factor_with_tube)[0]

    data_older_without_cone = make_data_by_run("mineral oil background with and without cone.txt", -90, 90, intensity_factor_older)[0]
    data_older_with_cone = make_data_by_run("mineral oil background with and without cone.txt", -90, 90, intensity_factor_older)[1]

    points_old_walls = make_points(data_old_walls[0], 0, 90, 1.47399, 0.5, data_old_walls[1], 405,
                                   photodiode_solid_angle, photodiode_angular_width, "cardboard walls")
    points_new_walls = make_points(data_new_walls[0], 0, 90, 1.47399, 0.5, data_new_walls[1], 405,
                                   photodiode_solid_angle, photodiode_angular_width, "new walls")
    points_with_tube = make_points(data_with_tube[0], 0, 90, 1.47399, 0.5, data_with_tube[1], 405,
                                   photodiode_solid_angle, photodiode_angular_width, "new walls with cone of darkness")
    points_older_without_cone = make_points(data_older_without_cone[0], 0, 90, 1.47399, 0.5, data_older_without_cone[1],
                                            405, photodiode_solid_angle, photodiode_angular_width, "older without cone")
    points_older_with_cone = make_points(data_older_with_cone[0], 0, 90, 1.47399, 0.5, data_older_with_cone[1],
                                            405, photodiode_solid_angle, photodiode_angular_width, "older with cone")


all_points = points_old_walls + points_new_walls + points_with_tube


#plot_points(points_old_walls + points_new_walls + points_with_tube, "Newer Background Tests", log=False, show=False)
#plot_points(points_older_with_cone + points_older_without_cone, "Older Background Tests", log=False, show=False)
#plot_points(points_older_with_cone + points_older_without_cone + points_new_walls + points_with_tube, "Old Cone vs New Cone", log=False, show=False)
plot_points(points_older_without_cone + points_new_walls, "Old vs New Without Cone", log=False, show=False)
plot_points(points_older_with_cone + points_with_tube, "Old vs New With Cone", log=False, show=False)
plt.show()
