# 3/6/2018

import numpy as np
import matplotlib.pyplot as plt
from plotting import make_data_by_run, plot_with_restarts

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
    flux_i = 0.00618 * 100e-6

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

    full_data = make_data_by_run("bubble tests all 3_12_18.txt", -90, 90, intensity_factor)

    full_data.reverse()

test_1 = True
if test_1:
    test_1_data = full_data[:]
    before_bubbles = test_1_data[:4]
    after_bubbles = test_1_data[4:6]
    after_sonic = test_1_data[6:]

    before_x = []
    before_y = []
    for data in before_bubbles:
        x_data = data[0]
        for x in x_data:
            before_x.append(x)
        y_data = data[1]
        for y in y_data:
            before_y.append(y)

    after_x = []
    after_y = []
    for data in after_bubbles:
        x_data = data[0]
        for x in x_data:
            after_x.append(x)
        y_data = data[1]
        for y in y_data:
            after_y.append(y)
			
    sonic_x = []
    sonic_y = []
    for data in after_sonic:
        x_data = data[0]
        for x in x_data:
            sonic_x.append(x)
        y_data = data[1]
        for y in y_data:
            sonic_y.append(y)

    plt.figure()
    plot_with_restarts(before_x, before_y, "before bubbles", "b")
    plot_with_restarts(after_x, after_y, "after bubbles", "r")
    plot_with_restarts(sonic_x, sonic_y, "after sonication", "g")

    #plt.yscale("log")
    plt.legend(loc=2)
    plt.xlabel("viewing angle (degrees)")
    plt.ylabel("intensity (flux/str)/(input flux)")
    #plt.title("Test 1")

plt.show()

