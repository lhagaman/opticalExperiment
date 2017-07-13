import numpy as np
import matplotlib.pyplot as plt
from Point import Point
import TSTR_fit


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

# solid angle of photodiode:
# 5.435 inches from sample
# 9 mm diameter
# in inches
distance_from_sample_to_photodiode = 5.435

# mm / (mm / inch)
photodiode_radius = (9 / 2.0) / 25.4

photodiode_solid_angle = np.pi * np.power(photodiode_radius, 2) / np.power(distance_from_sample_to_photodiode, 2)


make_air_points = True
if make_air_points:

    # volts * amps/volt
    flux_i = 0.00186 * 100e-6

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

    dir = "cylindrical_cell_no_polarizer_6-28/"

    data_30 = np.loadtxt(dir + "30 deg.txt", skiprows=1)
    data_45 = np.loadtxt(dir + "45 deg.txt", skiprows=1)
    data_55 = np.loadtxt(dir + "55 deg.txt", skiprows=1)

    data_30_x = [np.round(point[0], 2) for point in data_30]
    data_30_y = [intensity_factor * point[1] for point in data_30]

    data_45_x = [np.round(point[0], 2) for point in data_45]
    data_45_y = [intensity_factor * point[1] for point in data_45]

    data_55_x = [np.round(point[0], 2) for point in data_55]
    data_55_y = [intensity_factor * point[1] for point in data_55]

    points_30 = make_points(data_30_x, 0, 30, 1, 0.5, data_30_y, 405, photodiode_solid_angle,
                                     "30 degrees")
    points_45 = make_points(data_45_x, 0, 45, 1, 0.5, data_45_y, 405, photodiode_solid_angle,
                                     "45 degrees")
    points_55 = make_points(data_55_x, 0, 60, 1, 0.5, data_55_y, 405, photodiode_solid_angle,
                                     "55 degrees")

    points = points_30 + points_45 + points_55

    # cutting off regions which are dominated by reflections other than sample
    # off glass for low angles, and by background reflection for angles past parallel to sample surface
    cutoff_points_30 = []
    cutoff_data_30_x = []
    cutoff_data_30_y = []
    for point in points_30:
        if 10 < point.theta_r_in_degrees < 80:
            cutoff_points_30.append(point)
            cutoff_data_30_x.append(point.theta_r_in_degrees)
            cutoff_data_30_y.append(point.intensity)

    cutoff_points_45 = []
    cutoff_data_45_x = []
    cutoff_data_45_y = []
    for point in points_45:
        if 0 < point.theta_r_in_degrees < 80:
            cutoff_points_45.append(point)
            cutoff_data_45_x.append(point.theta_r_in_degrees)
            cutoff_data_45_y.append(point.intensity)

    cutoff_points_55 = []
    cutoff_data_55_x = []
    cutoff_data_55_y = []
    for point in points_55:
        if -20 < point.theta_r_in_degrees < 80:
            cutoff_points_55.append(point)
            cutoff_data_55_x.append(point.theta_r_in_degrees)
            cutoff_data_55_y.append(point.intensity)

    all_cutoff_points = cutoff_points_30 + cutoff_points_45 + cutoff_points_55

make_water_points = True
if make_water_points:

    # volts * amps/volt
    # estimated, changed for different angles of mirror in tube
    flux_i = 0.004 * 50e-6

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

    dir = "water_6-28/"

    w_data_30 = np.loadtxt(dir + "30 deg.txt", skiprows=1)
    w_data_45 = np.loadtxt(dir + "45 deg.txt", skiprows=1)
    w_data_55 = np.loadtxt(dir + "55 deg.txt", skiprows=1)

    w_data_30_x = [np.round(point[0], 2) for point in w_data_30]
    w_data_30_y = [intensity_factor * point[1] for point in w_data_30]

    w_data_45_x = [np.round(point[0], 2) for point in w_data_45]
    w_data_45_y = [intensity_factor * point[1] for point in w_data_45]

    w_data_55_x = [np.round(point[0], 2) for point in w_data_55]
    w_data_55_y = [intensity_factor * point[1] for point in w_data_55]

    plot_w_data = False

    if plot_w_data:
        plt.plot(w_data_30_x, w_data_30_y, label="30 degrees")
        plt.plot(w_data_45_x, w_data_45_y, label="45 degrees")
        plt.plot(w_data_55_x, w_data_55_y, label="55 degrees")

        plt.title("water, no polarizer")
        plt.legend()

        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.xlabel("viewing angle (degrees)")
        plt.show()

    w_points_30 = make_points(w_data_30_x, 0, 30, 1.33, 0.5, w_data_30_y, 405, photodiode_solid_angle,
                            "30 degrees")
    w_points_45 = make_points(w_data_45_x, 0, 45, 1.33, 0.5, w_data_45_y, 405, photodiode_solid_angle,
                            "45 degrees")
    w_points_55 = make_points(w_data_55_x, 0, 55, 1.33, 0.5, w_data_55_y, 405, photodiode_solid_angle,
                            "55 degrees")

    w_points = w_points_30 + w_points_45 + w_points_55

    # cutting off regions which are dominated by reflections other than sample
    # off glass for low angles, and by background reflection for angles past parallel to sample surface
    w_cutoff_points_30 = []
    w_cutoff_data_30_x = []
    w_cutoff_data_30_y = []
    for point in w_points_30:
        if 10 < point.theta_r_in_degrees < 80:
            w_cutoff_points_30.append(point)
            w_cutoff_data_30_x.append(point.theta_r_in_degrees)
            w_cutoff_data_30_y.append(point.intensity)

    w_cutoff_points_45 = []
    w_cutoff_data_45_x = []
    w_cutoff_data_45_y = []
    for point in w_points_45:
        if 0 < point.theta_r_in_degrees < 80:
            w_cutoff_points_45.append(point)
            w_cutoff_data_45_x.append(point.theta_r_in_degrees)
            w_cutoff_data_45_y.append(point.intensity)

    w_cutoff_points_55 = []
    w_cutoff_data_55_x = []
    w_cutoff_data_55_y = []
    for point in w_points_55:
        if -20 < point.theta_r_in_degrees < 80:
            w_cutoff_points_55.append(point)
            w_cutoff_data_55_x.append(point.theta_r_in_degrees)
            w_cutoff_data_55_y.append(point.intensity)

    all_w_cutoff_points = w_cutoff_points_30 + w_cutoff_points_45 + w_cutoff_points_55

make_mineral_oil_points = True
if make_mineral_oil_points:

    # volts * amps/volt
    # estimated, changed for different angles of mirror in tube
    flux_i = 0.0045 * 50e-6

    # amps per volt during measurements
    sensitivity = 100 * 1e-9

    # product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
    # intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
    intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

    dir = "mineral_oil_6-28/"

    m_data_30 = np.loadtxt(dir + "30 deg.txt", skiprows=1)
    m_data_45 = np.loadtxt(dir + "45 deg.txt", skiprows=1)
    m_data_55 = np.loadtxt(dir + "55 deg.txt", skiprows=1)

    m_data_30_x = [np.round(point[0], 2) for point in m_data_30]
    m_data_30_y = [intensity_factor * point[1] for point in m_data_30]

    m_data_45_x = [np.round(point[0], 2) for point in m_data_45]
    m_data_45_y = [intensity_factor * point[1] for point in m_data_45]

    m_data_55_x = [np.round(point[0], 2) for point in m_data_55]
    m_data_55_y = [intensity_factor * point[1] for point in m_data_55]

    plot_m_data = False

    if plot_m_data:
        plt.plot(m_data_30_x, m_data_30_y, label="30 degrees")
        plt.plot(m_data_45_x, m_data_45_y, label="45 degrees")
        plt.plot(m_data_55_x, m_data_55_y, label="55 degrees")

        plt.title("mineral oil")
        plt.legend()

        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.xlabel("viewing angle (degrees)")
        plt.show()

    m_points_30 = make_points(m_data_30_x, 0, 30, 1.461, 0.5, m_data_30_y, 405, photodiode_solid_angle,
                            "30 degrees")
    m_points_45 = make_points(m_data_45_x, 0, 45, 1.461, 0.5, m_data_45_y, 405, photodiode_solid_angle,
                            "45 degrees")
    m_points_55 = make_points(m_data_55_x, 0, 55, 1.461, 0.5, m_data_55_y, 405, photodiode_solid_angle,
                            "55 degrees")

    m_points = m_points_30 + m_points_45 + m_points_55

    # cutting off regions which are dominated by reflections other than sample
    # off glass for low angles, and by background reflection for angles past parallel to sample surface
    m_cutoff_points_30 = []
    m_cutoff_data_30_x = []
    m_cutoff_data_30_y = []
    for point in m_points_30:
        if 10 < point.theta_r_in_degrees < 80:
            m_cutoff_points_30.append(point)
            m_cutoff_data_30_x.append(point.theta_r_in_degrees)
            m_cutoff_data_30_y.append(point.intensity)

    m_cutoff_points_45 = []
    m_cutoff_data_45_x = []
    m_cutoff_data_45_y = []
    for point in m_points_45:
        if 0 < point.theta_r_in_degrees < 80:
            m_cutoff_points_45.append(point)
            m_cutoff_data_45_x.append(point.theta_r_in_degrees)
            m_cutoff_data_45_y.append(point.intensity)

    m_cutoff_points_55 = []
    m_cutoff_data_55_x = []
    m_cutoff_data_55_y = []
    for point in m_points_55:
        if -20 < point.theta_r_in_degrees < 80:
            m_cutoff_points_55.append(point)
            m_cutoff_data_55_x.append(point.theta_r_in_degrees)
            m_cutoff_data_55_y.append(point.intensity)

    all_m_cutoff_points = m_cutoff_points_30 + m_cutoff_points_45 + m_cutoff_points_55

print('fitting parameters...')

TSTR_parameters_30 = TSTR_fit.fit_parameters_and_angle(cutoff_points_30)
TSTR_parameters_45 = TSTR_fit.fit_parameters_and_angle(cutoff_points_45)
TSTR_parameters_55 = TSTR_fit.fit_parameters_and_angle(cutoff_points_55)

w_TSTR_parameters_30 = TSTR_fit.fit_parameters_and_angle(w_cutoff_points_30)
w_TSTR_parameters_45 = TSTR_fit.fit_parameters_and_angle(w_cutoff_points_45)
w_TSTR_parameters_55 = TSTR_fit.fit_parameters_and_angle(w_cutoff_points_55)

m_TSTR_parameters_30 = TSTR_fit.fit_parameters_and_angle(m_cutoff_points_30)
m_TSTR_parameters_45 = TSTR_fit.fit_parameters_and_angle(m_cutoff_points_45)
m_TSTR_parameters_55 = TSTR_fit.fit_parameters_and_angle(m_cutoff_points_55)

plt.scatter(TSTR_parameters_30[0], TSTR_fit.reflectance(TSTR_parameters_30[0], 1, 0.5, TSTR_parameters_30[1:]), marker='o', color='r')
plt.scatter(TSTR_parameters_45[0], TSTR_fit.reflectance(TSTR_parameters_45[0], 1, 0.5, TSTR_parameters_45[1:]), marker='o', color='r')
plt.scatter(TSTR_parameters_55[0], TSTR_fit.reflectance(TSTR_parameters_55[0], 1, 0.5, TSTR_parameters_55[1:]), marker='o', color='r')

plt.scatter(w_TSTR_parameters_30[0], TSTR_fit.reflectance(w_TSTR_parameters_30[0], 1.33, 0.5, w_TSTR_parameters_30[1:]), marker='o', color='g')
plt.scatter(w_TSTR_parameters_45[0], TSTR_fit.reflectance(w_TSTR_parameters_45[0], 1.33, 0.5, w_TSTR_parameters_45[1:]), marker='o', color='g')
plt.scatter(w_TSTR_parameters_55[0], TSTR_fit.reflectance(w_TSTR_parameters_55[0], 1.33, 0.5, w_TSTR_parameters_55[1:]), marker='o', color='g')

plt.scatter(m_TSTR_parameters_30[0], TSTR_fit.reflectance(m_TSTR_parameters_30[0], 1.461, 0.5, m_TSTR_parameters_30[1:]), marker='o', color='b')
plt.scatter(m_TSTR_parameters_45[0], TSTR_fit.reflectance(m_TSTR_parameters_45[0], 1.461, 0.5, m_TSTR_parameters_45[1:]), marker='o', color='b')
plt.scatter(m_TSTR_parameters_55[0], TSTR_fit.reflectance(m_TSTR_parameters_55[0], 1.461, 0.5, m_TSTR_parameters_55[1:]), marker='o', color='b')

print('integrating total reflectances...')

reflectance_30_total = []
reflectance_45_total = []
reflectance_55_total = []

w_reflectance_30_total = []
w_reflectance_45_total = []
w_reflectance_55_total = []

m_reflectance_30_total = []
m_reflectance_45_total = []
m_reflectance_55_total = []

theta_i_plot_list = [10, 20, 30, 40, 50, 60, 70, 80]

for extra_theta_i in [TSTR_parameters_30[0], TSTR_parameters_45[0], TSTR_parameters_55[0],
                      w_TSTR_parameters_30[0], w_TSTR_parameters_45[0], w_TSTR_parameters_55[0],
                      m_TSTR_parameters_30[0], m_TSTR_parameters_45[0], m_TSTR_parameters_55[0]]:
    if extra_theta_i not in theta_i_plot_list:
        theta_i_plot_list.append(extra_theta_i)
theta_i_plot_list = sorted(theta_i_plot_list)

print("    Theta_i list: " + str(theta_i_plot_list))

for theta_i in theta_i_plot_list:
    print('    ' + str(theta_i) + ' degrees...')

    reflectance_30_total.append(TSTR_fit.reflectance(theta_i, 1, 0.5, TSTR_parameters_30[1:]))
    reflectance_45_total.append(TSTR_fit.reflectance(theta_i, 1, 0.5, TSTR_parameters_45[1:]))
    reflectance_55_total.append(TSTR_fit.reflectance(theta_i, 1, 0.5, TSTR_parameters_55[1:]))

    w_reflectance_30_total.append(TSTR_fit.reflectance(theta_i, 1.33, 0.5, w_TSTR_parameters_30[1:]))
    w_reflectance_45_total.append(TSTR_fit.reflectance(theta_i, 1.33, 0.5, w_TSTR_parameters_45[1:]))
    w_reflectance_55_total.append(TSTR_fit.reflectance(theta_i, 1.33, 0.5, w_TSTR_parameters_55[1:]))

    m_reflectance_30_total.append(TSTR_fit.reflectance(theta_i, 1.461, 0.5, m_TSTR_parameters_30[1:]))
    m_reflectance_45_total.append(TSTR_fit.reflectance(theta_i, 1.461, 0.5, m_TSTR_parameters_45[1:]))
    m_reflectance_55_total.append(TSTR_fit.reflectance(theta_i, 1.461, 0.5, m_TSTR_parameters_55[1:]))

print('    finished')

"""
parameter_string_30 = 'air, fitted_theta_i: ' + str(TSTR_parameters_30[0]) + ' rho_L: ' + str(TSTR_parameters_30[1]) + ', n: ' + str(TSTR_parameters_30[2]) + ', gamma: ' + str(TSTR_parameters_30[3])
parameter_string_45 = 'air, fitted_theta_i: ' + str(TSTR_parameters_45[0]) + ' rho_L: ' + str(TSTR_parameters_45[1]) + ', n: ' + str(TSTR_parameters_45[2]) + ', gamma: ' + str(TSTR_parameters_45[3])
parameter_string_55 = 'air, fitted_theta_i: ' + str(TSTR_parameters_55[0]) + ' rho_L: ' + str(TSTR_parameters_55[1]) + ', n: ' + str(TSTR_parameters_55[2]) + ', gamma: ' + str(TSTR_parameters_55[3])
"""

plt.plot(theta_i_plot_list, reflectance_30_total, label='air', color='r')
plt.plot(theta_i_plot_list, reflectance_45_total, color='r')
plt.plot(theta_i_plot_list, reflectance_55_total, color='r')

plt.plot(theta_i_plot_list, w_reflectance_30_total, label='water', color='g')
plt.plot(theta_i_plot_list, w_reflectance_45_total, color='g')
plt.plot(theta_i_plot_list, w_reflectance_55_total, color='g')

plt.plot(theta_i_plot_list, m_reflectance_30_total, label='mineral oil', color='b')
plt.plot(theta_i_plot_list, m_reflectance_45_total, color='b')
plt.plot(theta_i_plot_list, m_reflectance_55_total, color='b')

plt.legend()
plt.title('Predicted Total Reflectance vs. Incident Angle for different materials')
plt.xlabel('incident angle (degrees)')
plt.ylabel('reflectance')
plt.show()
