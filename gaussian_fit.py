# using fit from chapter 6 of Levy thesis
# assumes theta_tilt (phi_r with other convention) is zero

import numpy as np
import scipy.optimize


# takes an array, easy to plot
# returns [[specular1, diffuse1], [specular2, diffuse2], ...]
def BRIDF_plotter(theta_PMT_in_degrees_array, theta_rot_in_degrees, parameters):
    theta_rot = theta_rot_in_degrees * np.pi / 180
    return_array = []
    for theta_rot_in_degrees in theta_PMT_in_degrees_array:
        theta_PMT = np.pi * theta_rot_in_degrees / 180
        return_array.append(BRIDF(theta_rot, theta_PMT, parameters))
    return return_array


def BRIDF_specular_plotter(theta_PMT_in_degrees_array, theta_rot_in_degrees, parameters):
    theta_rot = theta_rot_in_degrees * np.pi / 180
    return_array = []
    for theta_rot_in_degrees in theta_PMT_in_degrees_array:
        theta_PMT = np.pi * theta_rot_in_degrees / 180
        return_array.append(BRIDF_specular(theta_rot, theta_PMT, parameters))
    return return_array


def BRIDF_diffuse_plotter(theta_PMT_in_degrees_array, theta_rot_in_degrees, parameters):
    theta_rot = theta_rot_in_degrees * np.pi / 180
    return_array = []
    for theta_rot_in_degrees in theta_PMT_in_degrees_array:
        theta_PMT = np.pi * theta_rot_in_degrees / 180
        return_array.append(BRIDF_diffuse(theta_rot, theta_PMT, parameters))
    return return_array


# uses radians
# returns [specular, diffuse]
# only works for phi_r = 0
# R_1 and R_2 are functions of theta_rot
def BRIDF_pair(theta_PMT, theta_rot, parameters):
    sigma = parameters[0]
    R_1 = np.abs(parameters[1])
    R_2 = np.abs(parameters[2])
    alpha = np.arccos(np.sin(theta_rot) * np.sin(theta_PMT) + np.cos(theta_rot) * np.cos(theta_PMT))
    return [R_1 * np.exp(-np.power(alpha, 2) / (2 * np.power(sigma, 2))), R_2 * np.cos(theta_PMT)]


def BRIDF(theta_PMT, theta_rot, parameters):
    pair = BRIDF_pair(theta_PMT, theta_rot, parameters)
    return pair[0] + pair[1]


# used for other models to take this model as a piece
def BRIDF_all_parameters(theta_r, phi_r, theta_i, n_0, polarization, photodiode_solid_angle, parameters):
    pair = BRIDF_pair(theta_r, theta_i, parameters)
    return pair[0] + pair[1]


def BRIDF_specular(theta_PMT, theta_rot, parameters):
    return BRIDF_pair(theta_PMT, theta_rot, parameters)[0]


def BRIDF_diffuse(theta_PMT, theta_rot, parameters):
    return BRIDF_pair(theta_PMT, theta_rot, parameters)[1]


def unvectorized_fitter(independent_variables, sigma, R_1, R_2):
    theta_PMT = independent_variables[0] * np.pi / 180
    theta_rot = independent_variables[1] * np.pi / 180
    parameters = [sigma, R_1, R_2]
    return BRIDF(theta_PMT, theta_rot, parameters)


# each independent_variables has the form [theta_PMT_in_degrees, theta_rot_in_degrees]
def fitter(independent_variables_array, sigma, R_1, R_2):
    arr = []
    for independent_variables in independent_variables_array:
        arr.append(unvectorized_fitter(independent_variables, sigma, R_1, R_2))
    return arr


# each point has form [theta_PMT_in_degrees, phi_r_in_degrees, theta_rot_in_degrees, n_0, polarization, intensity]
# we only care about theta_PMT_in_degrees, theta_rot_in_degrees, intensity
def fit_parameters(points):
    independent_variables_array = []
    intensity_array = []
    for point in points:
        independent_variables_array.append([point.theta_r_in_degrees, point.theta_i_in_degrees])
        intensity_array.append(point.intensity)
    # initial parameters are approximately the ones found in the paper
    fit_params = scipy.optimize.curve_fit(fitter, independent_variables_array, intensity_array,
                                          p0=[0.1, 100, 50])[0]
    return [fit_params[0], fit_params[1], fit_params[2]]


def change_theta_i(points, new_theta_i):
    new_points = points[:]
    for point in new_points:
        point.theta_i_in_degrees = new_theta_i
    return new_points


def fit_parameters_and_angle(points):
    intensities = [point.intensity for point in points]
    # determine if there is a sharp peak
    max_point = points[intensities.index(max(intensities))]
    theta_r_peak = max_point.theta_r_in_degrees
    # peak is within 15 degrees of where it's expected
    if np.abs(max_point.theta_r_in_degrees - max_point.theta_i_in_degrees) < 15:
        # use peak as theta_i
        points_1 = change_theta_i(points, theta_r_peak)
        parameters_with_theta_i_peak = fit_parameters(points_1)
        theta_r_in_degrees_array = [point.theta_r_in_degrees for point in points]
        point_0 = points_1[0]
        fit_array = BRIDF_plotter(theta_r_in_degrees_array, point_0.theta_i_in_degrees, parameters_with_theta_i_peak)
        fit_peak_index = fit_array.index(max(fit_array))
        fit_peak = points[fit_peak_index].theta_r_in_degrees
        peak_offset = fit_peak - theta_r_peak
        points_2 = change_theta_i(points, theta_r_peak - peak_offset)
        return [theta_r_peak - peak_offset] + fit_parameters(points_2)

    else:
        # use experimental as theta_i
        return [max_point.theta_i_in_degrees] + fit_parameters(points)
