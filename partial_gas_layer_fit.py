# model made from scratch assuming that there is a gas layer next to the PTFE
# assumes that gas layer is much thicker than wavelength of light
# can take any other model for the reflection inside the gas layer

# I think I'm abandoning this because the total internal reflection component seems to always be much, much more
# reflective than in experiment

import numpy as np
import scipy.optimize
import gaussian_fit


# polarized electric field perpendicular to plane of incidence
def F_s(theta_i, n_0, n):
    # total internal reflection
    if n_0 / n * np.sin(theta_i) > 1:
        return 1
    theta_t = np.arcsin(n_0 / n * np.sin(theta_i))
    return np.power((n_0 * np.cos(theta_i) - n * np.cos(theta_t)) /
                    (n_0 * np.cos(theta_i) + n * np.cos(theta_t)), 2)


# polarized electric field parallel to plane of incidence
def F_p(theta_i, n_0, n):
    # total internal reflection
    if n_0 / n * np.sin(theta_i) > 1:
        return 1
    theta_t = np.arcsin(n_0 / n * np.sin(theta_i))
    return np.power((n_0 * np.cos(theta_t) - n * np.cos(theta_i)) /
                    (n_0 * np.cos(theta_t) + n * np.cos(theta_i)), 2)


def F_unpolarized(theta_i, n_0, n):
    return 0.5 * (F_s(theta_i, n_0, n) + F_p(theta_i, n_0, n))


def F(theta_i, n_0, n, polarization):
    return polarization * F_p(theta_i, n_0, n) + (1 - polarization) * F_s(theta_i, n_0, n)


# x is fraction of gas layer coverage (0 means no gas, 1 means gas layer everywhere)
def BRIDF(BRIDF_gas, theta_r, phi_r, theta_i, n_gas, n_liquid, polarization, photodiode_solid_angle, parameters, x):
    if not (0 <= x <= 1):
        return 1000000

    return x * BRIDF_with_gas_layer(BRIDF_gas, theta_r, phi_r, theta_i, n_gas, n_liquid, polarization,
                                    photodiode_solid_angle, parameters) + \
           (1 - x) * BRIDF_gas(theta_r, phi_r, theta_i, n_gas, polarization, photodiode_solid_angle,
                               parameters)


# BRIDF_gas is a total BRIDF function that we assume applies inside the gas layer
# BRIDF_gas takes all independent variables even if model only uses some of those
# assumes all reflected light has only one reflection off PTFE (no extended bouncing between edges of gas layer)
def BRIDF_with_gas_layer(BRIDF_gas, theta_r, phi_r, theta_i, n_gas, n_liquid, polarization, photodiode_solid_angle, parameters):

    # from liquid/gas interface reflection
    photodiode_angle = np.sqrt(photodiode_solid_angle)
    if np.abs(theta_r - theta_i) <= photodiode_angle and np.abs(phi_r) <= photodiode_angle:
        specular_component = F(theta_i, n_liquid, n_gas, polarization) / photodiode_solid_angle
    else:
        specular_component = 0

    # total internal reflection within liquid
    if n_liquid / n_gas * np.sin(theta_i) > 1:
        print(n_liquid / n_gas * np.sin(theta_i))
        print("theta_r: ", theta_r * 180 / np.pi, "not getting in gas layer")
        print("should be one: ", F(theta_i, n_liquid, n_gas, polarization))
        return specular_component

    theta_i_gas = np.arcsin(n_liquid / n_gas * np.sin(theta_i))

    # total internal reflection within gas, only reflection from liquid/gas interface contributes
    # should never happen because n_gas < n_liquid
    if n_gas / n_liquid * np.sin(theta_r) > 1:
        print("this shouldn't print")
        return specular_component

    theta_r_gas = np.arcsin(n_gas / n_liquid * np.sin(theta_r))

    frac_transmitted_into_gas = 1 - F(theta_i, n_liquid, n_gas, polarization)
    frac_transmitted_from_gas = 1 - F(theta_r_gas, n_liquid, n_gas, polarization)

    return specular_component + frac_transmitted_into_gas * frac_transmitted_from_gas * \
        BRIDF_gas(theta_r_gas, phi_r, theta_i_gas, n_gas, polarization, photodiode_solid_angle, parameters)


def BRIDF_plotter(BRIDF_gas, theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_gas, n_liquid, polarization,
                  photodiode_solid_angle, parameters, x):
    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        theta_r = np.pi * theta_r_in_degrees / 180
        return_array.append(BRIDF(BRIDF_gas, theta_r, phi_r, theta_i, n_gas, n_liquid, polarization,
                                  photodiode_solid_angle, parameters, x))
    return return_array


# independent variables has the form
# [theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, photodiode_solid_angle]
def unvectorized_fitter_gaussian(independent_variables, sigma, R_1, R_2, x):
    theta_r = independent_variables[0] * np.pi / 180
    phi_r = independent_variables[1] * np.pi / 180
    theta_i = independent_variables[2] * np.pi / 180

    n_gas = 1.3
    # from Measurement of the Refractive Index and Attenuation Length
    # of Liquid Xenon for its Scintillation Light, Solonov et. al.
    n_liquid = 1.69
    polarization = independent_variables[4]
    photodiode_solid_angle = independent_variables[5]
    parameters = [sigma, R_1, R_2]
    return BRIDF(gaussian_fit.BRIDF_all_parameters, theta_r, phi_r, theta_i, n_gas, n_liquid, polarization,
                 photodiode_solid_angle, parameters, x)


# takes arrays of independent variable lists
def fitter_gaussian(independent_variables_array, sigma, R_1, R_2, x):
    arr = []
    for independent_variables in independent_variables_array:
        arr.append(unvectorized_fitter_gaussian(independent_variables, sigma, R_1, R_2, x))
    return arr


# each point has form
# [theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity]
def fit_parameters_gaussian(points):
    independent_variables_array = []
    intensity_array = []
    for point in points:
        independent_variables_array.append([point.theta_r_in_degrees, point.phi_r_in_degrees, point.theta_i_in_degrees,
                                            point.n_0, point.polarization, point.photodiode_solid_angle])
        intensity_array.append(point.intensity)
    # initial parameters are the ones found in the paper
    fit_params = scipy.optimize.curve_fit(fitter_gaussian, independent_variables_array, intensity_array,
                                        p0=[0.1, 100, 50, 0.5])[0]
    return [fit_params[0], fit_params[1], fit_params[2], fit_params[3]]
