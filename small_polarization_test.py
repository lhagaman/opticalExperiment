# vertical no sample: 0.000829939 * 100e-6 amps

# horizontal no sample: 0.002304 * 100e-6 amps


import numpy as np
import matplotlib.pyplot as plt
from Point import Point
import TSTR_fit
from plotting import plot_with_TSTR_fit, make_points, make_data_by_run, plot_with_TSTR_fit_and_fitted_angles, plot_points

# in inches
distance_from_sample_to_photodiode = 5.435

# these are guessed, need to measure later
photodiode_height = 0.450
photodiode_width = 0.075

photodiode_angular_width_radians = 0.075 / 5.435
photodiode_angular_width = photodiode_angular_width_radians * 180. / np.pi

photodiode_area = photodiode_height * photodiode_width

photodiode_solid_angle = photodiode_area / np.power(distance_from_sample_to_photodiode, 2)

make_all_points = True
if make_all_points:
    # these are all in air through the tube, they correct for the polarizer and possible polarization of the laser
    flux_i_with_polarizer_horizontal = 0.002304 * 100e-6
    flux_i_no_polarizer = 0.00205 * 100e-6
    polarizer_transmission_horizontal = flux_i_with_polarizer_horizontal / flux_i_no_polarizer

    flux_i_with_polarizer_vertical = 0.000829939 * 100e-6
    flux_i_no_polarizer = 0.00205 * 100e-6
    polarizer_transmission_vertical = flux_i_with_polarizer_vertical / flux_i_no_polarizer

    # volts * amps/volt
    flux_i = 0.00205 * 100e-6
    # amps per volt during measurements
    sensitivity = 100 * 1e-9
    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor_air_vertical = sensitivity / (photodiode_solid_angle * flux_i * polarizer_transmission_vertical)
    intensity_factor_air_horizontal = sensitivity / (photodiode_solid_angle * flux_i * polarizer_transmission_horizontal)

    # volts * amps/volt
    flux_i = 0.00212 * 100e-6
    # amps per volt during measurements
    sensitivity = 100 * 1e-9
    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor_water_vertical = polarizer_transmission_vertical * sensitivity / (photodiode_solid_angle * flux_i)
    intensity_factor_water_horizontal = polarizer_transmission_vertical * sensitivity / (photodiode_solid_angle * flux_i)

    # volts * amps/volt
    flux_i = 0.00238 * 100e-6
    # amps per volt during measurements
    sensitivity = 100 * 1e-9
    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor_mineral_oil_vertical = polarizer_transmission_vertical * sensitivity / (photodiode_solid_angle * flux_i)
    intensity_factor_mineral_oil_horizontal = polarizer_transmission_vertical * sensitivity / (photodiode_solid_angle * flux_i)

    data_horizontal = make_data_by_run("45_deg_polarization_test.txt", 0, 80, intensity_factor_air_horizontal)[0]
    data_vertical = make_data_by_run("45_deg_polarization_test.txt", 0, 80, intensity_factor_air_vertical)[1]

    points_horizontal = make_points(data_horizontal[0], 0, 45, 1, 1, data_horizontal[1], 405, photodiode_solid_angle, "horizontal")
    points_vertical = make_points(data_vertical[0], 0, 45, 1, 0, data_vertical[1], 405, photodiode_solid_angle, "vertical")

plot = True
if plot:
    plot_with_TSTR_fit(points_vertical + points_horizontal, "Different Polarizations in Air, 45 degrees Without Slit")

