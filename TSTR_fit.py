# TSTR_fit

# specular follows a torrance-sparrow model
# diffuse follows a trowbridge-reitz distribution of surface normals combined with lambertian diffusion

# described in Reflectance of Polytetrafluoroethylene (PTFE) for Xenon_Scintillation Light, LIP-Coimbra

# called geometrical optical approximation (GOA) model in other Coimbra paper,
# but this doesn't include the correction factor N

import numpy as np
import scipy.optimize
import copy

# following equations could be refined for mu !~ mu_0, but teflon has mu~mu_0


# this function is used in calculating the total reflectance
def G_calc(theta_r, phi_r, theta_i, n_0, polarization, parameters):
    rho_L = parameters[0]
    n = parameters[1]
    gamma = parameters[2]

    theta_i_prime = 0.5 * np.arccos(np.cos(theta_i) * np.cos(theta_r) -
                                    np.sin(theta_i) * np.sin(theta_r) * np.cos(phi_r))

    theta_r_prime = 0.5 * np.arccos(np.cos(theta_i) * np.cos(theta_r) -
                                    np.sin(theta_i) * np.sin(theta_r))

    def G_prime(theta):
        return 2 / (1 + np.sqrt(1 + np.power(gamma * np.tan(theta), 2)))

        # this has negative of the inside of the H functions from the paper, I think it was a typo

    G = H(np.pi / 2 - theta_i_prime) * H(np.pi / 2 - theta_r_prime) * \
        G_prime(theta_i) * G_prime(theta_r)

    return G


# polarized electric field perpendicular to plane of incidence (vertical)
def F_s(theta_i, n_0, n):
    if np.abs(n_0 / n * np.sin(theta_i)) > 1:
        # total internal reflection
        return 1
    theta_t = np.arcsin(n_0 / n * np.sin(theta_i))
    return np.power((n_0 * np.cos(theta_i) - n * np.cos(theta_t)) /
                    (n_0 * np.cos(theta_i) + n * np.cos(theta_t)), 2)


# polarized electric field parallel to plane of incidence (horizontal)
def F_p(theta_i, n_0, n):
    if np.abs(n_0 / n * np.sin(theta_i)) > 1:
        # total internal reflection
        return 1
    theta_t = np.arcsin(n_0 / n * np.sin(theta_i))
    return np.power((n_0 * np.cos(theta_t) - n * np.cos(theta_i)) /
                    (n_0 * np.cos(theta_t) + n * np.cos(theta_i)), 2)


def F_unpolarized(theta_i, n_0, n):
    return 0.5 * (F_s(theta_i, n_0, n) + F_p(theta_i, n_0, n))


def F(theta_i, n_0, n, polarization):
    return polarization * F_p(theta_i, n_0, n) + (1 - polarization) * F_s(theta_i, n_0, n)


# heaviside step function
def H(x):
    return 0.5 * (np.sign(x) + 1)


# BRIDF as described in paper
# takes an array, easy to plot
# parameters has form [rho_L, n, gamma]


def BRIDF_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, parameters, average_angle=0):
    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        theta_r = np.pi * theta_r_in_degrees / 180
        return_array.append(BRIDF(theta_r, phi_r, theta_i, n_0, polarization, parameters))
    if average_angle == 0:
        return return_array
    else:
        return average_by_angle(theta_r_in_degrees_array, return_array, average_angle)


def average_by_angle(theta_r_in_degrees_array, intensity_array, angle):
    grouped = [[] for i in intensity_array]
    for i in range(len(intensity_array)):
        for j in range(len(intensity_array)):
            if np.abs(theta_r_in_degrees_array[j] - theta_r_in_degrees_array[i]) < angle / 2.:
                grouped[i].append(intensity_array[j])
    averaged = [np.average(grouped[i]) for i in range(len(intensity_array))]
    return averaged

def BRIDF_radiance_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, parameters):
    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        theta_r = np.pi * theta_r_in_degrees / 180
        return_array.append(BRIDF(theta_r, phi_r, theta_i, n_0, polarization, parameters) / np.cos(theta_r))
    return return_array


def BRIDF_specular_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees,
                           n_0, polarization, parameters):
    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        theta_r = np.pi * theta_r_in_degrees / 180
        return_array.append(BRIDF_specular(theta_r, phi_r, theta_i, n_0, polarization, parameters))
    return return_array


def BRIDF_diffuse_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees,
                          n_0, polarization, parameters):
    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        theta_r = np.pi * theta_r_in_degrees / 180
        return_array.append(BRIDF_diffuse(theta_r, phi_r, theta_i, n_0, polarization, parameters))
    return return_array


def BRIDF_pair(theta_r, phi_r, theta_i, n_0, polarization, parameters):

    rho_L = parameters[0]
    n = parameters[1]
    gamma = parameters[2]

    # from page 85 of coimbra thesis
    theta_prime = 0.5 * np.arccos(np.cos(theta_i) * np.cos(theta_r) -
        np.sin(theta_i) * np.sin(theta_r) * np.cos(phi_r))

    theta_i_prime = theta_prime
    theta_r_prime = theta_prime

    def P(alpha_):
        return np.power(gamma, 2) / \
               (np.pi * np.power(np.cos(alpha_), 4) *
                np.power(np.power(gamma, 2) + np.power(np.tan(alpha_), 2), 2))

    if theta_i == theta_r:
        alpha_specular = 0
    else:
        alpha_specular = np.arccos((np.cos(theta_i) + np.cos(theta_r)) / (2 * np.cos(theta_i_prime)))

    P_ = P(alpha_specular)

    if np.abs(n_0 / n * np.sin(theta_r)) > 1:
        # the W expression would have an arcsin error
        W = 0
    else:
        # There was a typo here (?), I changed by moving parentheses and taking a reciprocal
        W = (1 - F(theta_i, n_0, n, polarization)) * (1 - F_unpolarized(np.arcsin(n_0 / n * np.sin(theta_r)), n, n_0))

    A = 0.5 * np.power(gamma, 2) / (np.power(gamma, 2) + 0.92)

    theta_m = min(theta_i, theta_r)

    theta_M = max(theta_i, theta_r)

    B = 0.45 * np.power(gamma, 2) / (np.power(gamma, 2) + 0.25) * \
        H(np.cos(phi_r)) * np.cos(phi_r) * np.sin(theta_M) * np.tan(theta_m)

    def G_prime(theta):
        return 2 / (1 + np.sqrt(1 + np.power(gamma * np.tan(theta), 2)))

    # this has negative of the inside of the H functions from the paper, I think it was a typo
    G = H(np.pi / 2 - theta_i_prime) * H(np.pi / 2 - theta_r_prime) * \
        G_prime(theta_i) * G_prime(theta_r)

    F_ = F(theta_i_prime, n_0, n, polarization)

    return [F_ * G * P_ / (4 * np.cos(theta_i)),
        rho_L / np.pi * W * (1 - A + B) * np.cos(theta_r)]


def BRIDF_specular(theta_r, phi_r, theta_i, n_0, polarization, parameters):
    return BRIDF_pair(theta_r, phi_r, theta_i, n_0, polarization, parameters)[0]


def BRIDF_diffuse(theta_r, phi_r, theta_i, n_0, polarization, parameters):
    return BRIDF_pair(theta_r, phi_r, theta_i, n_0, polarization, parameters)[1]


def BRIDF(theta_r, phi_r, theta_i, n_0, polarization, parameters):
    pair = BRIDF_pair(theta_r, phi_r, theta_i, n_0, polarization, parameters)
    return pair[0] + pair[1]


# independent variables has the form [theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization]
def unvectorized_fitter(independent_variables, log_rho_L, log_n_minus_one, log_gamma):
    theta_r = independent_variables[0] * np.pi / 180
    phi_r = independent_variables[1] * np.pi / 180
    theta_i = independent_variables[2] * np.pi / 180
    n_0 = independent_variables[3]
    polarization = independent_variables[4]
    rho_L = np.exp(log_rho_L)
    n = np.exp(log_n_minus_one) + 1
    gamma = np.exp(log_gamma)
    parameters = [rho_L, n, gamma]
    return BRIDF(theta_r, phi_r, theta_i, n_0, polarization, parameters)


# independent variables without angle has the form [theta_r_in_degrees, phi_r_in_degrees, n_0, polarization]
def unvectorized_fitter_with_angle(independent_variables_without_angle, theta_i, log_rho_L, log_n_minus_one, log_gamma):
    theta_r = independent_variables_without_angle[0] * np.pi / 180
    phi_r = independent_variables_without_angle[1] * np.pi / 180
    n_0 = independent_variables_without_angle[2]
    polarization = independent_variables_without_angle[2]
    rho_L = np.exp(log_rho_L)
    n = np.exp(log_n_minus_one) + 1
    gamma = np.exp(log_gamma)
    parameters = [rho_L, n, gamma]
    return BRIDF(theta_r, phi_r, theta_i, n_0, polarization, parameters)


# takes arrays of independent variable lists
def fitter(independent_variables_array, log_rho_L, log_n_minus_one, log_gamma):
    arr = []
    average_angle = independent_variables_array[0][-1]
    for independent_variables in independent_variables_array:
        arr.append(unvectorized_fitter(independent_variables[:-1], log_rho_L, log_n_minus_one, log_gamma))
    if average_angle == 0:
        return arr
    else:
        theta_r_in_degrees_array = [ind_vars[0] for ind_vars in independent_variables_array]
	    # need to have theta_r in degrees! ind. vars is a list, each entry of which contains a list including theta_r (first)
        return average_by_angle(theta_r_in_degrees_array, arr, average_angle)

def fitter_with_angle(independent_variables_without_angle_array, theta_i, log_rho_L, log_n_minus_one, log_gamma):
    arr = []
    for independent_variables_without_angle in independent_variables_without_angle_array:
        fit_val = unvectorized_fitter_with_angle(independent_variables_without_angle, theta_i, log_rho_L, log_n_minus_one, log_gamma)
        arr.append(fit_val)
    return arr


def fit_parameters(points, p0=[0.5, 1.5, 0.05], average_angle=0):
    p0_log = [np.log(p0[0]), np.log(p0[1]-1), np.log(p0[2])]
    independent_variables_array = []
    intensity_array = []
    for point in points:
        independent_variables_array.append([point.theta_r_in_degrees, point.phi_r_in_degrees,
                                       point.theta_i_in_degrees, point.n_0, point.polarization, average_angle])
        intensity_array.append(point.intensity)
    # initial parameters are the ones found in the paper
    # fitter_p0 = fitter(independent_variables_array, p0_log[0], p0_log[1], p0_log[2])
    # res_fit_p0 = np.sum((np.array(intensity_array) - np.array(fitter_p0))**2)
    fit_params = scipy.optimize.curve_fit(fitter, independent_variables_array, intensity_array,
                                          p0=p0_log)[0]
    # fitter_popt = fitter(independent_variables_array, fit_params[0], fit_params[1], fit_params[2])
    # res_fit_popt = np.sum((np.array(intensity_array) - np.array(fitter_popt))**2)
    # print(res_fit_p0, res_fit_popt)
    return [np.exp(fit_params[0]), np.exp(fit_params)[1] + 1, np.exp(fit_params[2])]


# fits parameters, uses standard deviations
def fit_std(points, average_angle=0):
    independent_variables_array = []
    intensity_array = []
    std_array = []
    for point in points:
        independent_variables_array.append([point.theta_r_in_degrees, point.phi_r_in_degrees, point.theta_i_in_degrees, point.n_0, point.polarization, average_angle])
        intensity_array.append(point.intensity)
        std_array.append(point.std)
    # initial parameters are the ones found in the paper
    fit_params = scipy.optimize.curve_fit(fitter, independent_variables_array, intensity_array, p0=[np.log(0.5), np.log(1.5 - 1), np.log(0.05)], sigma=std_array)[0]
    return [np.exp(fit_params[0]), np.exp(fit_params)[1] + 1, np.exp(fit_params[2])]


def change_theta_i(points, new_theta_i):
    new_points = copy.deepcopy(points)
    for point in new_points:
        point.theta_i_in_degrees = new_theta_i
    return new_points


# returns [fitted_theta_i, rho_L, n, gamma]
# not intended for running data with more than one incident angle at a time
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
        fit_array = BRIDF_plotter(theta_r_in_degrees_array, point_0.phi_r_in_degrees, point_0.theta_i_in_degrees,
                                     point_0.n_0, point_0.polarization, parameters_with_theta_i_peak)
        fit_peak_index = fit_array.index(max(fit_array))
        fit_peak = points[fit_peak_index].theta_r_in_degrees
        peak_offset = fit_peak - theta_r_peak
        points_2 = change_theta_i(points, theta_r_peak - peak_offset)
        return [theta_r_peak - peak_offset] + fit_parameters(points_2)

    else:
        # use experimental as theta_i
        return [max_point.theta_i_in_degrees] + fit_parameters(points)


# returns [fitted_theta_i, rho_L, n, gamma]
# not intended for running data with more than one incident angle at a time
def fit_parameters_and_angle_with_starting_theta_i(points, starting_theta_i, average_angle=0):
    independent_variables_without_angle_array = []
    intensity_array = []
    for point in points:
        independent_variables_without_angle_array.append([point.theta_r_in_degrees, point.phi_r_in_degrees,
                                            point.n_0, point.polarization, average_angle])
        intensity_array.append(point.intensity)
    # initial parameters are the ones found in the paper
    fit_params = scipy.optimize.curve_fit(fitter_with_angle, independent_variables_without_angle_array, intensity_array,
                                          p0=[starting_theta_i, np.log(0.5), np.log(1.5 - 1), np.log(0.05)])[0]
    return [fit_params[0], np.exp(fit_params[1]), np.exp(fit_params)[2] + 1, np.exp(fit_params[3])]


def reflectance(theta_i_in_degrees, n_0, polarization, parameters):

    theta_i = theta_i_in_degrees * np.pi / 180

    a = lambda x: 0
    b = lambda x: np.pi / 2

    def BRIDF_int(theta_r, phi_r):
        # solid angle not used here
        G = G_calc(theta_r, phi_r, theta_i, n_0, polarization, parameters)
        return BRIDF(theta_r, phi_r, theta_i, n_0, polarization, parameters) * np.sin(theta_r) / G

    return scipy.integrate.dblquad(BRIDF_int, -np.pi, np.pi, a, b)[0]


def reflectance_diffuse(theta_i_in_degrees, n_0, polarization, parameters):

    theta_i = theta_i_in_degrees * np.pi / 180

    a = lambda x: 0
    b = lambda x: np.pi / 2

    def BRIDF_int(theta_r, phi_r):
        # solid angle not used here
        G = G_calc(theta_r, phi_r, theta_i, n_0, polarization, parameters)
        return BRIDF_diffuse(theta_r, phi_r, theta_i, n_0, polarization, parameters) * np.sin(theta_r) / G

    return scipy.integrate.dblquad(BRIDF_int, -np.pi, np.pi, a, b)[0]


def reflectance_specular(theta_i_in_degrees, n_0, polarization, parameters):

    theta_i = theta_i_in_degrees * np.pi / 180

    a = lambda x: 0
    b = lambda x: np.pi / 2

    def BRIDF_int(theta_r, phi_r):
        G = G_calc(theta_r, phi_r, theta_i, n_0, polarization, parameters)
        return BRIDF_specular(theta_r, phi_r, theta_i, n_0, polarization, parameters) * np.sin(theta_r) / G

    return scipy.integrate.dblquad(BRIDF_int, -np.pi, np.pi, a, b)[0]



