
# air through tube, 2.05 mV, 100uA/V
# water through tube, 2.12 mV, 100 uA/V
# mineral oil through tube: 2.38 mV, 100uA/V


import numpy as np
import matplotlib.pyplot as plt
from Point import Point
import TSTR_fit
from plotting import plot_TSTR_fit_one_set_of_parameters, plot_with_TSTR_fit, make_points, make_data_by_run, plot_with_TSTR_fit_and_fitted_angles, plot_points, plot_with_semi_empirical_fit

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

    flux_i_air_vertical = 0.0005400 * 100e-6
    flux_i_air_horizontal = 0.002304 * 100e-6

    flux_i_water_vertical = 0.0005135 * 100e-6
    flux_i_water_horizontal = 0.0017000 * 100e-6

    flux_i_mineral_oil_vertical = 0.000556 * 100e-6
    flux_i_mineral_oil_horizontal = 0.002330 * 100e-6

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor_air_vertical = sensitivity / (photodiode_solid_angle * flux_i_air_vertical)
    intensity_factor_air_horizontal = sensitivity / (photodiode_solid_angle * flux_i_air_horizontal)

    intensity_factor_water_vertical = sensitivity / (photodiode_solid_angle * flux_i_water_vertical)
    intensity_factor_water_horizontal = sensitivity / (photodiode_solid_angle * flux_i_water_horizontal)

    intensity_factor_mineral_oil_vertical = sensitivity / (photodiode_solid_angle * flux_i_mineral_oil_vertical)
    intensity_factor_mineral_oil_horizontal = sensitivity / (photodiode_solid_angle * flux_i_mineral_oil_horizontal)

    horizontal = True
    if horizontal:

        data_30_horizontal_air = make_data_by_run("new_polarization_comparison_data/horizontal/75, 60, 45, and 30 in air.txt", 0, 80, intensity_factor_air_horizontal)[0]
        data_45_horizontal_air = make_data_by_run("new_polarization_comparison_data/horizontal/75, 60, 45, and 30 in air.txt", 0, 80, intensity_factor_air_horizontal)[1]
        data_60_horizontal_air = make_data_by_run("new_polarization_comparison_data/horizontal/75, 60, 45, and 30 in air.txt", -20, 80, intensity_factor_air_horizontal)[2]
        data_75_horizontal_air = make_data_by_run("new_polarization_comparison_data/horizontal/75, 60, 45, and 30 in air.txt", -20, 80, intensity_factor_air_horizontal)[3]

        points_30_horizontal_air = make_points(data_30_horizontal_air[0], 0, 30, 1, 1, data_30_horizontal_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "30 degrees horizontal air")
        points_45_horizontal_air = make_points(data_45_horizontal_air[0], 0, 45, 1, 1, data_45_horizontal_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "45 degrees horizontal air")
        points_60_horizontal_air = make_points(data_60_horizontal_air[0], 0, 60, 1, 1, data_60_horizontal_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "60 degrees horizontal air")
        points_75_horizontal_air = make_points(data_75_horizontal_air[0], 0, 75, 1, 1, data_75_horizontal_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "75 degrees horizontal air")


        data_30_horizontal_water = make_data_by_run("new_polarization_comparison_data/horizontal/75, 60, 45, and 30 in water.txt", 0, 80, intensity_factor_water_horizontal)[0]
        data_45_horizontal_water = make_data_by_run("new_polarization_comparison_data/horizontal/75, 60, 45, and 30 in water.txt", 0, 80, intensity_factor_water_horizontal)[1]
        data_60_horizontal_water = make_data_by_run("new_polarization_comparison_data/horizontal/75, 60, 45, and 30 in water.txt", -20, 80, intensity_factor_water_horizontal)[2]
        data_75_horizontal_water = make_data_by_run("new_polarization_comparison_data/horizontal/75, 60, 45, and 30 in water.txt", -40, 80, intensity_factor_water_horizontal)[3]

        points_30_horizontal_water = make_points(data_30_horizontal_water[0], 0, 30, 1.33, 1, data_30_horizontal_water[1], 405, photodiode_solid_angle, photodiode_angular_width, "30 degrees horizontal water")
        points_45_horizontal_water = make_points(data_45_horizontal_water[0], 0, 45, 1.33, 1, data_45_horizontal_water[1], 405, photodiode_solid_angle, photodiode_angular_width, "45 degrees horizontal water")
        points_60_horizontal_water = make_points(data_60_horizontal_water[0], 0, 60, 1.33, 1, data_60_horizontal_water[1], 405, photodiode_solid_angle, photodiode_angular_width, "60 degrees horizontal water")
        points_75_horizontal_water = make_points(data_75_horizontal_water[0], 0, 75, 1.33, 1, data_75_horizontal_water[1], 405, photodiode_solid_angle, photodiode_angular_width, "75 degrees horizontal water")


        data_30_horizontal_mineral_oil = make_data_by_run("new_polarization_comparison_data/horizontal/75, 60, 45, and 30 in mineral oil.txt", 0, 80, intensity_factor_mineral_oil_horizontal)[0]
        data_45_horizontal_mineral_oil = make_data_by_run("new_polarization_comparison_data/horizontal/75, 60, 45, and 30 in mineral oil.txt", 0, 80, intensity_factor_mineral_oil_horizontal)[1]
        data_60_horizontal_mineral_oil = make_data_by_run("new_polarization_comparison_data/horizontal/75, 60, 45, and 30 in mineral oil.txt", -20, 80, intensity_factor_mineral_oil_horizontal)[2]
        data_75_horizontal_mineral_oil = make_data_by_run("new_polarization_comparison_data/horizontal/75, 60, 45, and 30 in mineral oil.txt", -40, 80, intensity_factor_mineral_oil_horizontal)[3]

        points_30_horizontal_mineral_oil = make_points(data_30_horizontal_mineral_oil[0], 0, 30, 1.461, 1, data_30_horizontal_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width, "30 degrees horizontal mineral oil")
        points_45_horizontal_mineral_oil = make_points(data_45_horizontal_mineral_oil[0], 0, 45, 1.461, 1, data_45_horizontal_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width, "45 degrees horizontal mineral oil")
        points_60_horizontal_mineral_oil = make_points(data_60_horizontal_mineral_oil[0], 0, 60, 1.461, 1, data_60_horizontal_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width, "60 degrees horizontal mineral oil")
        points_75_horizontal_mineral_oil = make_points(data_75_horizontal_mineral_oil[0], 0, 75, 1.461, 1, data_75_horizontal_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width, "75 degrees horizontal mineral oil")

    vertical = True
    if vertical:

        data_30_vertical_air = make_data_by_run("new_polarization_comparison_data/vertical/75, 60, 45, and 30 in air.txt", 0, 80, intensity_factor_air_vertical)[0]
        data_45_vertical_air = make_data_by_run("new_polarization_comparison_data/vertical/75, 60, 45, and 30 in air.txt", 0, 80, intensity_factor_air_vertical)[1]
        data_60_vertical_air = make_data_by_run("new_polarization_comparison_data/vertical/75, 60, 45, and 30 in air.txt", -20, 80, intensity_factor_air_vertical)[2]
        data_75_vertical_air = make_data_by_run("new_polarization_comparison_data/vertical/75, 60, 45, and 30 in air.txt", -20, 80, intensity_factor_air_vertical)[3]
        
        points_30_vertical_air = make_points(data_30_vertical_air[0], 0, 30, 1, 0, data_30_vertical_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "30 degrees vertical air")
        points_45_vertical_air = make_points(data_45_vertical_air[0], 0, 45, 1, 0, data_45_vertical_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "45 degrees vertical air")
        points_60_vertical_air = make_points(data_60_vertical_air[0], 0, 60, 1, 0, data_60_vertical_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "60 degrees vertical air")
        points_75_vertical_air = make_points(data_75_vertical_air[0], 0, 75, 1, 0, data_75_vertical_air[1], 405, photodiode_solid_angle, photodiode_angular_width, "75 degrees vertical air")
        

        data_30_vertical_water = make_data_by_run("new_polarization_comparison_data/vertical/75, 60, 45, and 30 in water.txt", 0, 80, intensity_factor_water_vertical)[0]
        data_45_vertical_water = make_data_by_run("new_polarization_comparison_data/vertical/75, 60, 45, and 30 in water.txt", 0, 80, intensity_factor_water_vertical)[1]
        data_60_vertical_water = make_data_by_run("new_polarization_comparison_data/vertical/75, 60, 45, and 30 in water.txt", -20, 80, intensity_factor_water_vertical)[2]
        data_75_vertical_water = make_data_by_run("new_polarization_comparison_data/vertical/75, 60, 45, and 30 in water.txt", -40, 80, intensity_factor_water_vertical)[3]

        points_30_vertical_water = make_points(data_30_vertical_water[0], 0, 30, 1.33, 0, data_30_vertical_water[1], 405, photodiode_solid_angle, photodiode_angular_width, "30 degrees vertical water")
        points_45_vertical_water = make_points(data_45_vertical_water[0], 0, 45, 1.33, 0, data_45_vertical_water[1], 405, photodiode_solid_angle, photodiode_angular_width, "45 degrees vertical water")
        points_60_vertical_water = make_points(data_60_vertical_water[0], 0, 60, 1.33, 0, data_60_vertical_water[1], 405, photodiode_solid_angle, photodiode_angular_width, "60 degrees vertical water")
        points_75_vertical_water = make_points(data_75_vertical_water[0], 0, 75, 1.33, 0, data_75_vertical_water[1], 405, photodiode_solid_angle, photodiode_angular_width, "75 degrees vertical water")


        data_30_vertical_mineral_oil = make_data_by_run("new_polarization_comparison_data/vertical/75, 60, 45, and 30 in mineral oil.txt", 0, 80, intensity_factor_mineral_oil_vertical)[0]
        data_45_vertical_mineral_oil = make_data_by_run("new_polarization_comparison_data/vertical/75, 60, 45, and 30 in water.txt", 0, 80, intensity_factor_mineral_oil_vertical)[1]
        data_60_vertical_mineral_oil = make_data_by_run("new_polarization_comparison_data/vertical/75, 60, 45, and 30 in water.txt", -20, 80, intensity_factor_mineral_oil_vertical)[2]
        data_75_vertical_mineral_oil = make_data_by_run("new_polarization_comparison_data/vertical/75, 60, 45, and 30 in water.txt", -40, 80, intensity_factor_mineral_oil_vertical)[3]

        points_30_vertical_mineral_oil = make_points(data_30_vertical_mineral_oil[0], 0, 30, 1.461, 0, data_30_vertical_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width, "30 degrees vertical mineral oil")
        points_45_vertical_mineral_oil = make_points(data_45_vertical_mineral_oil[0], 0, 45, 1.461, 0, data_45_vertical_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width, "45 degrees vertical mineral oil")
        points_60_vertical_mineral_oil = make_points(data_60_vertical_mineral_oil[0], 0, 60, 1.461, 0, data_60_vertical_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width, "60 degrees vertical mineral oil")
        points_75_vertical_mineral_oil = make_points(data_75_vertical_mineral_oil[0], 0, 75, 1.461, 0, data_75_vertical_mineral_oil[1], 405, photodiode_solid_angle, photodiode_angular_width, "75 degrees vertical mineral oil")


    points_air = points_30_horizontal_air + points_30_vertical_air + \
                 points_45_horizontal_air + points_45_vertical_air + \
                 points_60_horizontal_air + points_60_vertical_air + \
                 points_75_horizontal_air + points_75_vertical_air

    points_water = points_30_horizontal_water + points_30_vertical_water + \
                 points_45_horizontal_water + points_45_vertical_water + \
                 points_60_horizontal_water + points_60_vertical_water + \
                 points_75_horizontal_water + points_75_vertical_water

    points_mineral_oil = points_30_horizontal_mineral_oil + points_30_vertical_mineral_oil + \
                 points_45_horizontal_mineral_oil + points_45_vertical_mineral_oil + \
                 points_60_horizontal_mineral_oil + points_60_vertical_mineral_oil + \
                 points_75_horizontal_mineral_oil + points_75_vertical_mineral_oil

    all_points = points_air + points_water + points_mineral_oil

plot = True
if plot:

    fit = TSTR_fit.fit_parameters(points_30_horizontal_air + points_30_vertical_air + \
                 points_45_horizontal_air + points_45_vertical_air + \
                 points_60_horizontal_air + points_60_vertical_air)

    plot_TSTR_fit_one_set_of_parameters(points_30_horizontal_air + points_45_horizontal_air + points_60_horizontal_air, fit, "Air Horizontal")
    plot_TSTR_fit_one_set_of_parameters(points_30_vertical_air + points_45_vertical_air + points_60_vertical_air, fit, "Air Vertical")

    plot_TSTR_fit_one_set_of_parameters(points_30_horizontal_water + points_45_horizontal_water + points_60_horizontal_water, fit, "Water Horizontal")
    plot_TSTR_fit_one_set_of_parameters(points_30_vertical_water + points_45_vertical_water + points_60_vertical_water, fit, "Water Vertical")

    plot_TSTR_fit_one_set_of_parameters(points_30_horizontal_mineral_oil + points_45_horizontal_mineral_oil + points_60_horizontal_mineral_oil, fit, "Mineral Oil Horizontal")
    plot_TSTR_fit_one_set_of_parameters(points_30_vertical_mineral_oil + points_45_vertical_mineral_oil + points_60_vertical_mineral_oil, fit, "Mineral Oil Vertical")

    plt.show()

    #plot_with_TSTR_fit(points_30_vertical_air + points_30_horizontal_air, "Different Polarizations in Air With Slit")
    #plot_with_TSTR_fit(points_45_vertical_air + points_45_horizontal_air, "Different Polarizations in Air With Slit")
    #plot_with_TSTR_fit(points_45_vertical_air + points_45_horizontal_air, "Different Polarizations in Air With Slit")


    #plot_with_TSTR_fit(points_60_vertical_air + points_60_horizontal_air, "Different Polarizations in Air With Slit")
    #plot_with_TSTR_fit(points_75_vertical_air + points_75_horizontal_air, "Different Polarizations in Air With Slit")
    #plot_with_TSTR_fit(points_air, "Different Polarizations in Air With Slit")


    """
    #plot_points(points_30_vertical_water + points_30_horizontal_water, "Different Polarizations in water With Slit")
    #plot_points(points_45_vertical_water + points_45_horizontal_water, "Different Polarizations in water With Slit")
    #plot_with_TSTR_fit(points_30_vertical_water + points_30_horizontal_water, "Different Polarizations in water With Slit")
    #plot_with_TSTR_fit(points_45_vertical_water + points_45_horizontal_water, "Different Polarizations in water With Slit")
    plot_with_TSTR_fit(points_60_vertical_water + points_60_horizontal_water, "Different Polarizations in water With Slit")
    plot_with_TSTR_fit(points_75_vertical_water + points_75_horizontal_water, "Different Polarizations in water With Slit")
    plot_with_TSTR_fit(points_water, "Different Polarizations in water With Slit")
    """

    """
    plot_points(points_30_vertical_mineral_oil + points_30_horizontal_mineral_oil, "Different Polarizations in Mineral Oil With Slit")
    plot_points(points_45_vertical_mineral_oil + points_45_horizontal_mineral_oil, "Different Polarizations in Mineral Oil With Slit")
    plot_points(points_60_vertical_mineral_oil + points_60_horizontal_mineral_oil, "Different Polarizations in Mineral Oil With Slit")
    plot_points(points_75_vertical_mineral_oil + points_75_horizontal_mineral_oil, "Different Polarizations in Mineral Oil With Slit")
    plot_points(points_mineral_oil, "Different Polarizations in mineral_oil With Slit")
    """

    """

    plot_with_TSTR_fit(points_30_vertical_mineral_oil + points_30_horizontal_mineral_oil, "Different Polarizations in Mineral Oil With Slit")
    plot_with_TSTR_fit(points_45_vertical_mineral_oil + points_45_horizontal_mineral_oil, "Different Polarizations in Mineral Oil With Slit")
    plot_with_TSTR_fit(points_60_vertical_mineral_oil + points_60_horizontal_mineral_oil, "Different Polarizations in Mineral Oil With Slit")
    plot_with_TSTR_fit(points_75_vertical_mineral_oil + points_75_horizontal_mineral_oil, "Different Polarizations in Mineral Oil With Slit")
    plot_with_TSTR_fit(points_mineral_oil, "Different Polarizations in mineral_oil With Slit")
    """

    """
    plot_with_TSTR_fit(points_30_vertical_air, "Different Polarizations in air With Slit")
    plot_with_TSTR_fit(points_45_vertical_air, "Different Polarizations in air With Slit")
    plot_with_TSTR_fit(points_60_vertical_air, "Different Polarizations in air With Slit")
    plot_with_TSTR_fit(points_75_vertical_air, "Different Polarizations in air With Slit")

    plot_with_TSTR_fit(points_30_vertical_water, "Different Polarizations in Water With Slit")
    plot_with_TSTR_fit(points_45_vertical_water, "Different Polarizations in Water With Slit")
    plot_with_TSTR_fit(points_60_vertical_water, "Different Polarizations in Water With Slit")
    plot_with_TSTR_fit(points_75_vertical_water, "Different Polarizations in Water With Slit")

    plot_with_TSTR_fit(points_30_vertical_mineral_oil, "Different Polarizations in mineral_oil With Slit")
    plot_with_TSTR_fit(points_45_vertical_mineral_oil, "Different Polarizations in mineral_oil With Slit")
    plot_with_TSTR_fit(points_60_vertical_mineral_oil, "Different Polarizations in mineral_oil With Slit")
    plot_with_TSTR_fit(points_75_vertical_mineral_oil, "Different Polarizations in mineral_oil With Slit")
    """

    """
    plot_with_TSTR_fit(points_30_horizontal_air, "Different Polarizations in air With Slit")
    plot_with_TSTR_fit(points_45_horizontal_air, "Different Polarizations in air With Slit")
    plot_with_TSTR_fit(points_60_horizontal_air, "Different Polarizations in air With Slit")
    plot_with_TSTR_fit(points_75_horizontal_air, "Different Polarizations in air With Slit")

    plot_with_TSTR_fit(points_30_horizontal_water, "Different Polarizations in Water With Slit")
    plot_with_TSTR_fit(points_45_horizontal_water, "Different Polarizations in Water With Slit")
    plot_with_TSTR_fit(points_60_horizontal_water, "Different Polarizations in Water With Slit")
    plot_with_TSTR_fit(points_75_horizontal_water, "Different Polarizations in Water With Slit")

    plot_with_TSTR_fit(points_30_horizontal_mineral_oil, "Different Polarizations in mineral_oil With Slit")
    plot_with_TSTR_fit(points_45_horizontal_mineral_oil, "Different Polarizations in mineral_oil With Slit")
    plot_with_TSTR_fit(points_60_horizontal_mineral_oil, "Different Polarizations in mineral_oil With Slit")
    plot_with_TSTR_fit(points_75_horizontal_mineral_oil, "Different Polarizations in mineral_oil With Slit")
    """




