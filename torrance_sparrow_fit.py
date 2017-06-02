#TSTR_fit

# specular follows a torrance-sparrow model
# diffuse follows a trowbridge-reitz distribution of surface normals combined with lambertian diffusion

# detailed in Reflectance of Polytetrafluoroethylene (PTFE) for Xenon_Scintillation Light, LIP-Coimbra

import numpy as np
import scipy.optimize

# following equations could be refined for mu !~ mu_0, but teflon has mu~mu_0


# polarized electric field perpendicular to plane of incidence
def F_s(theta_i, n_0, n):
    theta_t = np.arcsin(n_0 / n * np.sin(theta_i))
    return np.power((n_0 * np.cos(theta_i) - n * np.cos(theta_t)) /
                    (n_0 * np.cos(theta_i) + n * np.cos(theta_t)), 2)


# polarized electric field parallel to plane of incidence
def F_p(theta_i, n_0, n):
    theta_t = np.arcsin(n_0 / n * np.sin(theta_i))
    return np.power((n_0 * np.cos(theta_t) - n * np.cos(theta_i)) /
                    (n_0 * np.cos(theta_t) + n * np.cos(theta_i)), 2)


def F_unpolarized(theta_i, n_0, n):
    return 0.5 * (F_s(theta_i, n_0, n) + F_p(theta_i, n_0, n))


# heaviside step function
def H(x):
    return 0.5 * (np.sign(x) + 1)


# BRIDF as described in paper
# takes an array, easy to plot
# parameters has form [rho_L, n, gamma]


def BRIDF_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, parameters):
    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        theta_r = np.pi * theta_r_in_degrees / 180
        return_array.append(BRIDF(theta_r, phi_r, theta_i, n_0, polarization, parameters))
    return return_array


def BRIDF_specular_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, parameters):
    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        theta_r = np.pi * theta_r_in_degrees / 180
        return_array.append(BRIDF_specular(theta_r, phi_r, theta_i, n_0, polarization, parameters))
    return return_array


def BRIDF_diffuse_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, parameters):
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
    theta_i_prime = 0.5 * np.arccos(np.cos(theta_i) * np.cos(theta_r) -
        np.sin(theta_i) * np.sin(theta_r) * np.cos(phi_r))

    # this part wasn't given in the paper, not entirely sure about it, just swapped i and r
    theta_r_prime = 0.5 * np.arccos(np.cos(theta_i) * np.cos(theta_r) -
        np.sin(theta_i) * np.sin(theta_r))

    def P(alpha_):
        return np.power(gamma, 2) / \
               (np.pi * np.power(np.cos(alpha_), 4) *
                np.power(np.power(gamma, 2) + np.power(np.tan(alpha_), 2), 2))

    alpha_specular = np.arccos((np.cos(theta_i) + np.cos(theta_r)) / (2 * np.cos(theta_i_prime)))

    P_ = P(alpha_specular)

    if polarization == 0:
        F = F_unpolarized
    elif polarization == 1:
        F = F_s
    elif polarization == 2:
        F = F_p
    else:
        print("INVALID POLARIZATION")

    # There was a typo here, I changed by moving parentheses and taking a reciprocal
    W = (1 - F(theta_i, n_0, n)) * (1 - F(np.arcsin(n_0 / n * np.sin(theta_r)), n_0, n))

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

    F_ = F(theta_i, n_0, n)

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


# takes arrays of independent variable lists
def fitter(independent_variables_array, log_rho_L, log_n_minus_one, log_gamma):
    arr = []
    for independent_variables in independent_variables_array:
        arr.append(unvectorized_fitter(independent_variables, log_rho_L, log_n_minus_one, log_gamma))
    return arr


# each point has form [theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity]
def fit_parameters(points):
    independent_variables_array = []
    intensity_array = []
    for point in points:
        independent_variables_array.append(point[:5])
        intensity_array.append(point[5])
    # initial parameters are the ones found in the paper
    fit_params = scipy.optimize.curve_fit(fitter, independent_variables_array, intensity_array,
                                          p0=[np.log(0.5), np.log(1.5 - 1), np.log(0.05)])[0]
    return [np.exp(fit_params[0]), np.exp(fit_params)[1] + 1, np.exp(fit_params[2])]
