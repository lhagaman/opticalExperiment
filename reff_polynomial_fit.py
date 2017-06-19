# this follows the paper Bi-directional Reflectance of Dry and Submerged Labsphere Spectralon Plaque
# by Hao Zhang and Kenneth J. Voss

# only suitable for incident light parallel to the normal to the sample
# uses degrees, not radians, throughout

import numpy as np
import scipy.optimize


def REFF(theta_r, parameters):
    c_1 = parameters[0]
    c_2 = parameters[1]
    c_3 = parameters[2]
    return c_1 + c_2 * theta_r * theta_r + c_3 * theta_r * theta_r * theta_r * theta_r


def BRIDF(theta_r, parameters):
    theta_i = 0
    mu_0 = np.cos(theta_i)
    return mu_0 / np.pi * REFF(theta_r, parameters)


def BRIDF_plotter(theta_r_array, parameters):
    return_array = []
    for theta_r in theta_r_array:
        return_array.append(BRIDF(theta_r, parameters))
    return return_array


def unvectorized_fitter(theta_r, c_1, c_2, c_3):
    parameters = [c_1, c_2, c_3]
    return BRIDF(theta_r, parameters)


def fitter(theta_r_array, c_1, c_2, c_3):
    arr = []
    for theta_r in theta_r_array:
        arr.append(unvectorized_fitter(theta_r, c_1, c_2, c_3))
    return arr


def fit_parameters(points):
    theta_r_array = []
    intensity_array = []
    for point in points:
        if point.theta_i_in_degrees != 0:
            print("not normal incidence, not expected to fit well")
        theta_r_array.append(point.theta_r_in_degrees)
        intensity_array.append(point.intensity)
    # initial parameters are approximately the ones found in the paper
    fit_params = scipy.optimize.curve_fit(fitter, theta_r_array, intensity_array,
                                          p0=[1, -1.52, -3.14])[0]
    return [fit_params[0], fit_params[1], fit_params[2]]
