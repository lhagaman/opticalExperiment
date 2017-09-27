
# air through tube, 2.05 mV, 100uA/V
# water through tube, 2.12 mV, 100 uA/V
# mineral oil through tube: 2.38 mV, 100uA/V


import numpy as np
import matplotlib.pyplot as plt
from Point import Point
import TSTR_fit
from plotting import subtract_background, plot_with_TSTR_fit, make_points, make_data_by_run, plot_with_TSTR_fit_and_fitted_angles, plot_points


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
    make_points_100 = True
    if make_points_100:
        # volts * amps/volt
        flux_i = 0.002582 * 100e-6

        # amps per volt during measurements
        sensitivity = 100 * 1e-9

        # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
        # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
        intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

        data_30 = make_data_by_run("glycerol_tests/glycerol 30,45,60,52.5,75.txt", -90, 90, intensity_factor)[4]
        data_45 = make_data_by_run("glycerol_tests/glycerol 30,45,60,52.5,75.txt", -90, 90, intensity_factor)[3]
        data_60 = make_data_by_run("glycerol_tests/glycerol 30,45,60,52.5,75.txt", -90, 90, intensity_factor)[2]
        data_52_5 = make_data_by_run("glycerol_tests/glycerol 30,45,60,52.5,75.txt", -90, 90, intensity_factor)[1]
        data_75 = make_data_by_run("glycerol_tests/glycerol 30,45,60,52.5,75.txt", -90, 90, intensity_factor)[0]

        points_100_30 = make_points(data_30[0], 0, 30, 1.47399, 0.5, data_30[1], 405, photodiode_solid_angle,
                                    photodiode_angular_width, "30 degrees 100% glycerol")
        points_100_45 = make_points(data_45[0], 0, 45, 1.47399, 0.5, data_45[1], 405, photodiode_solid_angle,
                                    photodiode_angular_width, "45 degrees 100% glycerol")
        points_100_52_5 = make_points(data_52_5[0], 0, 52.5, 1.47399, 0.5, data_52_5[1], 405, photodiode_solid_angle,
                                    photodiode_angular_width, "52.5 degrees 100% glycerol")
        points_100_60 = make_points(data_60[0], 0, 60, 1.47399, 0.5, data_60[1], 405, photodiode_solid_angle,
                                    photodiode_angular_width, "60 degrees 100% glycerol")
        points_100_75 = make_points(data_75[0], 0, 75, 1.47399, 0.5, data_75[1], 405, photodiode_solid_angle,
                                    photodiode_angular_width, "75 degrees 100% glycerol")

        points_100 = points_100_30 + points_100_45 + points_100_52_5 + points_100_60 + points_100_75

    make_points_0 = True
    if make_points_0:
        # volts * amps/volt
        flux_i = 0.002042 * 100e-6

        # amps per volt during measurements
        sensitivity = 100 * 1e-9

        # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
        # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
        intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

        data_30 = make_data_by_run("glycerol_tests/water 30,45,60,75.txt", -90, 90, intensity_factor)[3]
        data_45 = make_data_by_run("glycerol_tests/water 30,45,60,75.txt", -90, 90, intensity_factor)[2]
        data_60 = make_data_by_run("glycerol_tests/water 30,45,60,75.txt", -90, 90, intensity_factor)[1]
        data_75 = make_data_by_run("glycerol_tests/water 30,45,60,75.txt", -90, 90, intensity_factor)[0]

        points_0_30 = make_points(data_30[0], 0, 30, 1.33303, 0.5, data_30[1], 405, photodiode_solid_angle,
                                    photodiode_angular_width, "30 degrees 0% glycerol")
        points_0_45 = make_points(data_45[0], 0, 45, 1.33303, 0.5, data_45[1], 405, photodiode_solid_angle,
                                    photodiode_angular_width, "45 degrees 0% glycerol")
        points_0_60 = make_points(data_60[0], 0, 60, 1.33303, 0.5, data_60[1], 405, photodiode_solid_angle,
                                    photodiode_angular_width, "60 degrees 0% glycerol")
        points_0_75 = make_points(data_75[0], 0, 75, 1.33303, 0.5, data_75[1], 405, photodiode_solid_angle,
                                    photodiode_angular_width, "75 degrees 0% glycerol")

        points_0 = points_0_30 + points_0_45 + points_0_60 + points_0_75

    make_points_50 = True

    if make_points_50:
        # volts * amps/volt
        flux_i = 0.002515 * 100e-6

        # amps per volt during measurements
        sensitivity = 100 * 1e-9

        # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
        # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
        intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

        data_30 = make_data_by_run("glycerol_tests/50_ 30,45,60,52.5,75.txt", -90, 90, intensity_factor)[4]
        data_45 = make_data_by_run("glycerol_tests/50_ 30,45,60,52.5,75.txt", -90, 90, intensity_factor)[3]
        data_60 = make_data_by_run("glycerol_tests/50_ 30,45,60,52.5,75.txt", -90, 90, intensity_factor)[2]
        data_52_5 = make_data_by_run("glycerol_tests/50_ 30,45,60,52.5,75.txt", -90, 90, intensity_factor)[1]
        data_75 = make_data_by_run("glycerol_tests/50_ 30,45,60,52.5,75.txt", -90, 90, intensity_factor)[0]

        points_50_30 = make_points(data_30[0], 0, 30, 1.39809, 0.5, data_30[1], 405, photodiode_solid_angle,
                                  photodiode_angular_width, "30 degrees 50% glycerol")
        points_50_45 = make_points(data_45[0], 0, 45, 1.39809, 0.5, data_45[1], 405, photodiode_solid_angle,
                                  photodiode_angular_width, "45 degrees 50% glycerol")
        points_50_52_5 = make_points(data_52_5[0], 0, 52.5, 1.39809, 0.5, data_52_5[1], 405, photodiode_solid_angle,
                                    photodiode_angular_width, "52.5 degrees 50% glycerol")
        points_50_60 = make_points(data_60[0], 0, 60, 1.39809, 0.5, data_60[1], 405, photodiode_solid_angle,
                                  photodiode_angular_width, "60 degrees 50% glycerol")
        points_50_75 = make_points(data_75[0], 0, 75, 1.39809, 0.5, data_75[1], 405, photodiode_solid_angle,
                                  photodiode_angular_width, "75 degrees 50% glycerol")

        points_50 = points_50_30 + points_50_45 + points_50_52_5 + points_50_60 + points_50_75

    make_points_mineral_oil = True

    if make_points_mineral_oil:
        # volts * amps/volt
        flux_i = 0.002603 * 100e-6

        # amps per volt during measurements
        sensitivity = 100 * 1e-9

        # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
        # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
        intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

        data_30 = make_data_by_run("glycerol_tests/mineral oil 30,45,60,52.5,75.txt", -90, 90, intensity_factor)[4]
        data_45 = make_data_by_run("glycerol_tests/mineral oil 30,45,60,52.5,75.txt", -90, 90, intensity_factor)[3]
        data_60 = make_data_by_run("glycerol_tests/mineral oil 30,45,60,52.5,75.txt", -90, 90, intensity_factor)[2]
        data_52_5 = make_data_by_run("glycerol_tests/mineral oil 30,45,60,52.5,75.txt", -90, 90, intensity_factor)[1]
        data_75 = make_data_by_run("glycerol_tests/mineral oil 30,45,60,52.5,75.txt", -90, 90, intensity_factor)[0]

        points_mineral_oil_30 = make_points(data_30[0], 0, 30, 1.39809, 0.5, data_30[1], 405, photodiode_solid_angle,
                                            photodiode_angular_width, "30 degrees mineral oil")
        points_mineral_oil_45 = make_points(data_45[0], 0, 45, 1.39809, 0.5, data_45[1], 405, photodiode_solid_angle,
                                            photodiode_angular_width, "45 degrees mineral oil")
        points_mineral_oil_52_5 = make_points(data_52_5[0], 0, 52.5, 1.39809, 0.5, data_52_5[1], 405,
                                              photodiode_solid_angle,
                                              photodiode_angular_width, "52.5 degrees mineral oil")
        points_mineral_oil_60 = make_points(data_60[0], 0, 60, 1.39809, 0.5, data_60[1], 405, photodiode_solid_angle,
                                            photodiode_angular_width, "60 degrees mineral oil")
        points_mineral_oil_75 = make_points(data_75[0], 0, 75, 1.39809, 0.5, data_75[1], 405, photodiode_solid_angle,
                                            photodiode_angular_width, "75 degrees mineral oil")

        points_mineral_oil = points_mineral_oil_30 + points_mineral_oil_45 + points_mineral_oil_52_5 + points_mineral_oil_60 + points_mineral_oil_75

    make_background_points = True

    if make_background_points:
        # volts * amps/volt
        flux_i = 0.002603 * 100e-6

        # amps per volt during measurements
        sensitivity = 100 * 1e-9

        # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
        # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
        intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

        data = make_data_by_run("glycerol_tests/no sample background in water.txt", -90, 90, intensity_factor)[0]

        points_background = make_points(data[0], 0, 30, 1.39809, 0.5, data[1], 405, photodiode_solid_angle,
                                            photodiode_angular_width, "background data in mineral oil")

    make_mineral_oil_tape_points = True

    if make_mineral_oil_tape_points:
        # volts * amps/volt
        flux_i = 0.002695 * 100e-6

        # amps per volt during measurements
        sensitivity = 100 * 1e-9

        # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
        # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
        intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

        data_45 = make_data_by_run("glycerol_tests/mineral oil background with tape over sample 45,0,60,75.txt", -90, 90, intensity_factor)[3]
        data_0 = make_data_by_run("glycerol_tests/mineral oil background with tape over sample 45,0,60,75.txt", -90, 90, intensity_factor)[2]
        data_60 = make_data_by_run("glycerol_tests/mineral oil background with tape over sample 45,0,60,75.txt", -90, 90, intensity_factor)[1]
        data_75 = make_data_by_run("glycerol_tests/mineral oil background with tape over sample 45,0,60,75.txt", -90, 90, intensity_factor)[0]

        points_tape_45 = make_points(data_45[0], 0, 45, 1.39809, 0.5, data_45[1], 405, photodiode_solid_angle,
                                        photodiode_angular_width, "tape in mineral oil, 45 degrees")
        points_tape_0 = make_points(data_0[0], 0, 0, 1.39809, 0.5, data_0[1], 405, photodiode_solid_angle,
                                     photodiode_angular_width, "tape in mineral oil, 0 degrees")
        points_tape_60 = make_points(data_60[0], 0, 60, 1.39809, 0.5, data_60[1], 405, photodiode_solid_angle,
                                     photodiode_angular_width, "tape in mineral oil, 60 degrees")
        points_tape_75 = make_points(data_75[0], 0, 75, 1.39809, 0.5, data_75[1], 405, photodiode_solid_angle,
                                     photodiode_angular_width, "tape in mineral oil, 75 degrees")

        points_tape = points_tape_0 + points_tape_45 + points_tape_60 + points_tape_75

    make_9_26_background = True
    if make_9_26_background:
        # volts * amps/volt
        flux_i = 0.002638 * 100e-6

        # amps per volt during measurements
        sensitivity = 100 * 1e-9

        # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
        # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
        intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

        d_1 = make_data_by_run("glycerol_tests/9_26 backgrounds.txt", -90, 90, intensity_factor)[3]
        d_2 = make_data_by_run("glycerol_tests/9_26 backgrounds.txt", -90, 90, intensity_factor)[2]
        d_3 = make_data_by_run("glycerol_tests/9_26 backgrounds.txt", -90, 90, intensity_factor)[1]
        d_4 = make_data_by_run("glycerol_tests/9_26 backgrounds.txt", -90, 90, intensity_factor)[0]

        p_1 = make_points(d_1[0], 0, 90, 1.39809, 0.5, d_1[1], 405, photodiode_solid_angle,
                          photodiode_angular_width, "background in mineral oil with cell at usual angle")
        p_2 = make_points(d_2[0], 0, 90, 1.39809, 0.5, d_2[1], 405, photodiode_solid_angle,
                          photodiode_angular_width, "background in mineral oil with cell at different angle")
        p_3 = make_points(d_3[0], 0, 90, 1.39809, 0.5, d_3[1], 405, photodiode_solid_angle,
                          photodiode_angular_width, "background in mineral oil with cell at usual angle \nand black foam board angled upward behind cell")
        p_4 = make_points(d_4[0], 0, 90, 1.39809, 0.5, d_4[1], 405, photodiode_solid_angle,
                          photodiode_angular_width, "background in mineral oil with cell at usual angle \nand black foam board angled down, leaning against cell opposite photodiode")

        points_9_26 = p_1 + p_2 + p_3 + p_4

    make_cone_comparison = True
    if make_cone_comparison:
        # volts * amps/volt
        flux_i = 0.002603 * 100e-6

        # amps per volt during measurements
        sensitivity = 100 * 1e-9

        # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
        # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
        intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

        d_1 = make_data_by_run("glycerol_tests/mineral oil background with and without cone.txt", -90, 90, intensity_factor)[1]
        d_2 = make_data_by_run("glycerol_tests/mineral oil background with and without cone.txt", -90, 90, intensity_factor)[0]

        p_1 = make_points(d_1[0], 0, 90, 1.39809, 0.5, d_1[1], 405, photodiode_solid_angle,
                          photodiode_angular_width, "background in mineral oil with cone")
        p_2 = make_points(d_2[0], 0, 90, 1.39809, 0.5, d_2[1], 405, photodiode_solid_angle,
                          photodiode_angular_width, "background in mineral oil without cone")

        cone_comparison_points = p_1 + p_2


for p in points_tape:
    theta_rot = 90 - p.theta_r_in_degrees
    theta_r = 180 - p.theta_i_in_degrees - theta_rot
    p.theta_r_in_degrees = theta_r

#plot_points(points_tape, "Tape in Mineral Oil", log=False)

"""

points_mineral_oil_30_subtracted = subtract_background(points_mineral_oil_30, points_background)
points_mineral_oil_45_subtracted = subtract_background(points_mineral_oil_45, points_background)
points_mineral_oil_52_5_subtracted = subtract_background(points_mineral_oil_52_5, points_background)
points_mineral_oil_60_subtracted = subtract_background(points_mineral_oil_60, points_background)
points_mineral_oil_75_subtracted = subtract_background(points_mineral_oil_75, points_background)
"""


#plot_points(points_background, "", log=False, show=True)
"""
plot_points(points_mineral_oil_30 + points_mineral_oil_30_subtracted,
            "Mineral Oil 30 degrees and Background", log=False, show=False)
plot_points(points_mineral_oil_45 + points_mineral_oil_45_subtracted,
            "Mineral Oil 45 degrees and Background", log=False, show=False)
plot_points(points_mineral_oil_52_5 + points_mineral_oil_52_5_subtracted,
            "Mineral Oil 52.5 degrees and Background", log=False, show=False)
plot_points(points_mineral_oil_60 + points_mineral_oil_60_subtracted,
            "Mineral Oil 60 degrees and Background", log=False, show=False)
plot_points(points_mineral_oil_75 + points_mineral_oil_75_subtracted,
            "Mineral Oil 75 degrees and Background", log=False, show=False)

plt.show()

"""

#plot_points(points_9_26, "9/26 Background Measurements", log=False)

plot_points(cone_comparison_points, "Cone Of Darkness Comparison", log=False)


