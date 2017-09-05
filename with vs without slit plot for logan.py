# unknown normalization constants...

import numpy as np
import matplotlib.pyplot as plt
from Point import Point
import TSTR_fit
from plotting import plot_with_TSTR_fit, make_points, make_data_by_run, plot_with_TSTR_fit_and_fitted_angles, plot_points

# in inches
distance_from_sample_to_photodiode = 5.435

photodiode_height = 0.450
photodiode_width = 0.075

photodiode_angular_width_radians = 0.075 / 5.435
photodiode_angular_width = photodiode_angular_width_radians * 180. / np.pi

photodiode_area = photodiode_height * photodiode_width

photodiode_solid_angle = photodiode_area / np.power(distance_from_sample_to_photodiode, 2)

make_without_slit_points = True
if make_without_slit_points:
    photodiode_radius = (9 / 2.0) / 25.4

    photodiode_solid_angle = np.pi * np.power(photodiode_radius, 2) / np.power(distance_from_sample_to_photodiode, 2)

    photodiode_angular_width_radians = photodiode_radius / distance_from_sample_to_photodiode
    photodiode_angular_width = photodiode_angular_width_radians * 180. / np.pi

    photodiode_area = photodiode_height * photodiode_width

    # these are all in air through the tube, they correct for the polarizer and possible polarization of the laser

    # DIDN'T HAVE DATA FOR THESE, MADE UP IN ORDER TO COMPARE SHAPES
    flux_i_with_polarizer_horizontal_air = 0.0018 * 100e-6 * 5
    flux_i_with_polarizer_horizontal_water = 0.0018 * 100e-6 * 5
    flux_i_with_polarizer_horizontal_mineral_oil = 0.0018 * 100e-6 * 5

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i

    intensity_factor_horizontal_air = sensitivity / (photodiode_solid_angle * flux_i_with_polarizer_horizontal_air)
    intensity_factor_horizontal_water = sensitivity / (photodiode_solid_angle * flux_i_with_polarizer_horizontal_water)
    intensity_factor_horizontal_mineral_oil = sensitivity / (photodiode_solid_angle * flux_i_with_polarizer_horizontal_mineral_oil)

    horizontal = True
    if horizontal:
        data_30_horizontal_air = make_data_by_run("without_slit_7_18_horizontal/75, 60, 45, and 30 deg in air.txt", 0, 85, intensity_factor_horizontal_air)[0]
        data_45_horizontal_air = make_data_by_run("without_slit_7_18_horizontal/75, 60, 45, and 30 deg in air.txt", 0, 85, intensity_factor_horizontal_air)[1]
        data_60_horizontal_air = make_data_by_run("without_slit_7_18_horizontal/75, 60, 45, and 30 deg in air.txt", 0, 85, intensity_factor_horizontal_air)[2]
        data_75_horizontal_air = make_data_by_run("without_slit_7_18_horizontal/75, 60, 45, and 30 deg in air.txt", 0, 85, intensity_factor_horizontal_air)[3]

        points_30_horizontal_air = make_points(data_30_horizontal_air[0], 0, 30, 1, 1, data_30_horizontal_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "30 degrees horizontal air without slit")
        points_45_horizontal_air = make_points(data_45_horizontal_air[0], 0, 45, 1, 1, data_45_horizontal_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "45 degrees horizontal air without slit")
        points_60_horizontal_air = make_points(data_60_horizontal_air[0], 0, 60, 1, 1, data_60_horizontal_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "60 degrees horizontal air without slit")
        points_75_horizontal_air = make_points(data_75_horizontal_air[0], 0, 79, 1, 1, data_75_horizontal_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "75 degrees horizontal air without slit")


        data_30_horizontal_water = make_data_by_run("without_slit_7_18_horizontal/75, 60, 45, and 30 deg in water.txt", 0, 85, intensity_factor_horizontal_water)[0]
        data_45_horizontal_water = make_data_by_run("without_slit_7_18_horizontal/75, 60, 45, and 30 deg in water.txt", 0, 85, intensity_factor_horizontal_water)[1]
        data_60_horizontal_water = make_data_by_run("without_slit_7_18_horizontal/75, 60, 45, and 30 deg in water.txt", 0, 85, intensity_factor_horizontal_water)[2]
        data_75_horizontal_water = make_data_by_run("without_slit_7_18_horizontal/75, 60, 45, and 30 deg in water.txt", 0, 85, intensity_factor_horizontal_water)[3]

        points_30_horizontal_water = make_points(data_30_horizontal_water[0], 0, 30, 1.33, 1, data_30_horizontal_water[1], 405, photodiode_solid_angle, photodiode_angular_width, "30 degrees horizontal water without slit")
        points_45_horizontal_water = make_points(data_45_horizontal_water[0], 0, 45, 1.33, 1, data_45_horizontal_water[1], 405, photodiode_solid_angle, photodiode_angular_width, "45 degrees horizontal water without slit")
        points_60_horizontal_water = make_points(data_60_horizontal_water[0], 0, 60, 1.33, 1, data_60_horizontal_water[1], 405, photodiode_solid_angle, photodiode_angular_width, "60 degrees horizontal water without slit")
        points_75_horizontal_water = make_points(data_75_horizontal_water[0], 0, 75, 1.33, 1, data_75_horizontal_water[1], 405, photodiode_solid_angle, photodiode_angular_width, "75 degrees horizontal water without slit")


        data_30_horizontal_mineral_oil = make_data_by_run("without_slit_7_18_horizontal/75, 60, 45, and 30 deg in mineral oil.txt", 0, 85, intensity_factor_horizontal_mineral_oil)[0]
        data_45_horizontal_mineral_oil = make_data_by_run("without_slit_7_18_horizontal/75, 60, 45, and 30 deg in mineral oil.txt", 0, 85, intensity_factor_horizontal_mineral_oil)[1]
        data_60_horizontal_mineral_oil = make_data_by_run("without_slit_7_18_horizontal/75, 60, 45, and 30 deg in mineral oil.txt", 0, 85, intensity_factor_horizontal_mineral_oil)[2]
        data_75_horizontal_mineral_oil = make_data_by_run("without_slit_7_18_horizontal/75, 60, 45, and 30 deg in mineral oil.txt", 0, 85, intensity_factor_horizontal_mineral_oil)[3]

        points_30_horizontal_mineral_oil = make_points(data_30_horizontal_mineral_oil[0], 0, 30, 1.461, 1, data_30_horizontal_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width, "30 degrees horizontal mineral_oil without slit")
        points_45_horizontal_mineral_oil = make_points(data_45_horizontal_mineral_oil[0], 0, 45, 1.461, 1, data_45_horizontal_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width, "45 degrees horizontal mineral_oil without slit")
        points_60_horizontal_mineral_oil = make_points(data_60_horizontal_mineral_oil[0], 0, 60, 1.461, 1, data_60_horizontal_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width, "60 degrees horizontal mineral_oil without slit")
        points_75_horizontal_mineral_oil = make_points(data_75_horizontal_mineral_oil[0], 0, 75, 1.461, 1, data_75_horizontal_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width, "75 degrees horizontal mineral_oil without slit")

    points_air = points_30_horizontal_air + points_45_horizontal_air + points_60_horizontal_air + points_75_horizontal_air

    points_air_without_75 = points_30_horizontal_air + points_45_horizontal_air + points_60_horizontal_air

    points_water = points_30_horizontal_water + points_45_horizontal_water + points_60_horizontal_water + points_75_horizontal_water

    points_water_without_75 = points_30_horizontal_water + points_45_horizontal_water + points_60_horizontal_water

    points_mineral_oil = points_30_horizontal_mineral_oil + points_45_horizontal_mineral_oil + points_60_horizontal_mineral_oil + points_75_horizontal_mineral_oil

    points_mineral_oil_without_75 = points_30_horizontal_mineral_oil + points_45_horizontal_mineral_oil + points_60_horizontal_mineral_oil

    all_points = points_air + points_water + points_mineral_oil

make_with_slit_points = True
if make_with_slit_points:

    photodiode_height = 0.450
    photodiode_width = 0.075

    photodiode_angular_width_radians = 0.075 / 5.435
    photodiode_angular_width = photodiode_angular_width_radians * 180. / np.pi

    photodiode_area = photodiode_height * photodiode_width

    # these are all in air through the tube, they correct for the polarizer and possible polarization of the laser
    flux_i_with_polarizer_horizontal_air = 0.00205 * 100e-6
    flux_i_with_polarizer_horizontal_water = 0.00212 * 100e-6
    flux_i_with_polarizer_horizontal_mineral_oil = 0.00238 * 100e-6

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i

    intensity_factor_horizontal_air = sensitivity / (photodiode_solid_angle * flux_i_with_polarizer_horizontal_air)
    intensity_factor_horizontal_water = sensitivity / (photodiode_solid_angle * flux_i_with_polarizer_horizontal_water)
    intensity_factor_horizontal_mineral_oil = sensitivity / (photodiode_solid_angle * flux_i_with_polarizer_horizontal_mineral_oil)

    horizontal = True
    if horizontal:
        data_30_horizontal_air = make_data_by_run("horizontal_7_18/75, 60, 45, and 30 in air.txt", 0, 85, intensity_factor_horizontal_air)[0]
        data_45_horizontal_air = make_data_by_run("horizontal_7_18/75, 60, 45, and 30 in air.txt", 0, 85, intensity_factor_horizontal_air)[1]
        data_60_horizontal_air = make_data_by_run("horizontal_7_18/75, 60, 45, and 30 in air.txt", 0, 85, intensity_factor_horizontal_air)[2]
        data_75_horizontal_air = make_data_by_run("horizontal_7_18/75, 60, 45, and 30 in air.txt", 0, 85, intensity_factor_horizontal_air)[3]

        points_30_horizontal_air_slit = make_points(data_30_horizontal_air[0], 0, 30, 1, 1, data_30_horizontal_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "30 degrees horizontal air with slit")
        points_45_horizontal_air_slit = make_points(data_45_horizontal_air[0], 0, 45, 1, 1, data_45_horizontal_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "45 degrees horizontal air with slit")
        points_60_horizontal_air_slit = make_points(data_60_horizontal_air[0], 0, 60, 1, 1, data_60_horizontal_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "60 degrees horizontal air with slit")
        points_75_horizontal_air_slit = make_points(data_75_horizontal_air[0], 0, 79, 1, 1, data_75_horizontal_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "75 degrees horizontal air with slit")

        data_30_horizontal_water = make_data_by_run("horizontal_7_18/75, 60, 45, and 30 in water.txt", 0, 85, intensity_factor_horizontal_water)[0]
        data_45_horizontal_water = make_data_by_run("horizontal_7_18/75, 60, 45, and 30 in water.txt", 0, 85, intensity_factor_horizontal_water)[1]
        data_60_horizontal_water = make_data_by_run("horizontal_7_18/75, 60, 45, and 30 in water.txt", 0, 85, intensity_factor_horizontal_water)[2]
        data_75_horizontal_water = make_data_by_run("horizontal_7_18/75, 60, 45, and 30 in water.txt", 0, 85, intensity_factor_horizontal_water)[3]

        points_30_horizontal_water_slit = make_points(data_30_horizontal_water[0], 0, 30, 1.33, 1, data_30_horizontal_water[1], 405, photodiode_solid_angle, photodiode_angular_width,
                                                 "30 degrees horizontal water with slit")
        points_45_horizontal_water_slit = make_points(data_45_horizontal_water[0], 0, 45, 1.33, 1, data_45_horizontal_water[1], 405, photodiode_solid_angle, photodiode_angular_width,
                                                 "45 degrees horizontal water with slit")
        points_60_horizontal_water_slit = make_points(data_60_horizontal_water[0], 0, 60, 1.33, 1, data_60_horizontal_water[1], 405, photodiode_solid_angle, photodiode_angular_width,
                                                 "60 degrees horizontal water with slit")
        points_75_horizontal_water_slit = make_points(data_75_horizontal_water[0], 0, 75, 1.33, 1, data_75_horizontal_water[1], 405, photodiode_solid_angle, photodiode_angular_width,
                                                 "75 degrees horizontal water with slit")

        data_30_horizontal_mineral_oil = make_data_by_run("horizontal_7_18/75, 60, 45, and 30 in mineral oil.txt", 0, 85, intensity_factor_horizontal_mineral_oil)[0]
        data_45_horizontal_mineral_oil = make_data_by_run("horizontal_7_18/75, 60, 45, and 30 in mineral oil.txt", 0, 85, intensity_factor_horizontal_mineral_oil)[1]
        data_60_horizontal_mineral_oil = make_data_by_run("horizontal_7_18/75, 60, 45, and 30 in mineral oil.txt", 0, 85, intensity_factor_horizontal_mineral_oil)[2]
        data_75_horizontal_mineral_oil = make_data_by_run("horizontal_7_18/75, 60, 45, and 30 in mineral oil.txt", 0, 85, intensity_factor_horizontal_mineral_oil)[3]

        points_30_horizontal_mineral_oil_slit = make_points(data_30_horizontal_mineral_oil[0], 0, 30, 1.461, 1, data_30_horizontal_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width,
                                                       "30 degrees horizontal mineral_oil with slit")
        points_45_horizontal_mineral_oil_slit = make_points(data_45_horizontal_mineral_oil[0], 0, 45, 1.461, 1, data_45_horizontal_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width,
                                                       "45 degrees horizontal mineral_oil with slit")
        points_60_horizontal_mineral_oil_slit = make_points(data_60_horizontal_mineral_oil[0], 0, 60, 1.461, 1, data_60_horizontal_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width,
                                                       "60 degrees horizontal mineral_oil with slit")
        points_75_horizontal_mineral_oil_slit = make_points(data_75_horizontal_mineral_oil[0], 0, 75, 1.461, 1, data_75_horizontal_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width,
                                                       "75 degrees horizontal mineral_oil with slit")

    points_air_slit = points_30_horizontal_air_slit + points_45_horizontal_air_slit + points_60_horizontal_air_slit + points_75_horizontal_air_slit

    points_air_without_75_slit = points_30_horizontal_air_slit + points_45_horizontal_air_slit + points_60_horizontal_air_slit

    points_water_slit = points_30_horizontal_water_slit + points_45_horizontal_water_slit + points_60_horizontal_water_slit + points_75_horizontal_water_slit

    points_water_without_75_slit = points_30_horizontal_water_slit + points_45_horizontal_water_slit + points_60_horizontal_water_slit

    points_mineral_oil_slit = points_30_horizontal_mineral_oil_slit + points_45_horizontal_mineral_oil_slit + points_60_horizontal_mineral_oil_slit + points_75_horizontal_mineral_oil_slit

    points_mineral_oil_without_75_slit = points_30_horizontal_mineral_oil_slit + points_45_horizontal_mineral_oil_slit + points_60_horizontal_mineral_oil_slit

    all_points_slit = points_air_slit + points_water_slit + points_mineral_oil_slit

plot_air = False
if plot_air:
    #plot_points(points_air_without_75 + points_air_without_75_slit, "Horizontal Polarization In Air With And Without Slit")
    plot_points(points_30_horizontal_air + points_30_horizontal_air_slit, "Horizontal Polarization In Air With And Without Slit")
    plot_points(points_45_horizontal_air + points_45_horizontal_air_slit, "Horizontal Polarization In Air With And Without Slit")
    plot_points(points_60_horizontal_air + points_60_horizontal_air_slit, "Horizontal Polarization In Air With And Without Slit")
    plot_points(points_75_horizontal_air + points_75_horizontal_air_slit, "Horizontal Polarization In Air With And Without Slit")

plot_water = False
if plot_water:
    #plot_with_TSTR_fit(points_water, "Horizontal Polarization in Water Without Slit")
    plot_points(points_30_horizontal_water + points_30_horizontal_water_slit, "Horizontal Polarization In Water With And Without Slit")
    plot_points(points_45_horizontal_water + points_45_horizontal_water_slit, "Horizontal Polarization In Water With And Without Slit")
    plot_points(points_60_horizontal_water + points_60_horizontal_water_slit, "Horizontal Polarization In Water With And Without Slit")
    plot_points(points_75_horizontal_water + points_75_horizontal_water_slit, "Horizontal Polarization In Water With And Without Slit")

plot_mineral_oil = True
if plot_mineral_oil:
    #plot_with_TSTR_fit(points_mineral_oil_without_75, "Horizontal Polarization in Mineral Oil Without Slit")
    plot_points(points_30_horizontal_mineral_oil + points_30_horizontal_mineral_oil_slit, "Horizontal Polarization In Mineral Oil With And Without Slit")
    plot_points(points_45_horizontal_mineral_oil + points_45_horizontal_mineral_oil_slit, "Horizontal Polarization In Mineral Oil With And Without Slit")
    plot_points(points_60_horizontal_mineral_oil + points_60_horizontal_mineral_oil_slit, "Horizontal Polarization In Mineral Oil With And Without Slit")
    plot_points(points_75_horizontal_mineral_oil + points_75_horizontal_mineral_oil_slit, "Horizontal Polarization In Mineral Oil With And Without Slit")

