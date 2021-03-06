# Semi-Empirical model

# diffuse described by Wolff-Oren-Nayar model
# specular lobe described by geometrical optical approximation (GOA) model
# specular spike described by wave scattering Kirchhoff approximation

# described in Study of the Reflectance Distributions of Fluoropolymers
# and Other Rough Surfaces with Interest to Scintillation Detectors by Coimbra, 2009
# equation 4.109, page 121
# (called the "semi-empirical model" in paper)

import numpy as np
import scipy.optimize


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


def F(theta_i, n_0, n, polarization):
    return polarization * F_p(theta_i, n_0, n) + (1 - polarization) * F_s(theta_i, n_0, n)


# heaviside step function
def H(x):
    return 0.5 * (np.sign(x) + 1)


# we don't use the kappa parameter because mu~mu_0 in teflon
# parameters has the form [rho_L, n, K, gamma]
# returns [diffuse, specular lobe, specular spike]
def BRIDF_trio(theta_r, phi_r, theta_i, n_0, polarization, photodiode_angle, parameters):

    rho_L = parameters[0]
    n = parameters[1]
    K = parameters[2]
    gamma = parameters[3]

    theta_i_prime = 0.5 * np.arccos(np.cos(theta_i) * np.cos(theta_r) -
                                    np.sin(theta_i) * np.sin(theta_r) * np.cos(phi_r))

    theta_r_prime = 0.5 * np.arccos(np.cos(theta_i) * np.cos(theta_r) -
                                    np.sin(theta_i) * np.sin(theta_r))

    alpha_specular = np.arccos((np.cos(theta_i) + np.cos(theta_r)) / (2 * np.cos(theta_i_prime)))

    P = np.power(gamma, 2) / \
        (np.pi * np.power(np.cos(alpha_specular), 4) *
         np.power(np.power(gamma, 2) + np.power(np.tan(alpha_specular), 2), 2))

    F_ = F(theta_i, n_0, n, polarization)

    def G_prime(theta):
        return 2 / (1 + np.sqrt(1 + np.power(gamma * np.tan(theta), 2)))

    # this has negative of the inside of the H functions from the paper, I think it was a typo
    G = H(np.pi / 2 - theta_i_prime) * H(np.pi / 2 - theta_r_prime) * \
        G_prime(theta_i) * G_prime(theta_r)

    # assume that after scattering, light is unpolarized (page 115 of paper)
    W = (1 - F(theta_i, n_0, n, polarization)) * (1 - F_unpolarized(np.arcsin(n_0 / n * np.sin(theta_r)), n, n_0))

    # from page 152 of paper
    Lambda = np.exp(-K * np.cos(theta_i))

    gamma_squared = gamma * gamma

    # make the fitter dislike invalid gamma values
    if 1 - gamma_squared < 0:
        return [1000000, 1000000, 1000000]

    script_n_0 = 1 / (1 - gamma_squared) - \
                gamma_squared / (1 - gamma_squared) * \
                np.arctanh(np.sqrt(1 - gamma_squared)) / np.sqrt(1 - gamma_squared)

    script_n = gamma_squared / (2 * (1 - gamma_squared)) - \
                gamma_squared * (2 - gamma_squared) / (1 - gamma_squared) * \
                np.arctanh(np.sqrt(1 - gamma_squared)) / np.sqrt(1 - gamma_squared)

    N = G_prime(theta_i) * G_prime(theta_r) * \
        (script_n_0 - np.tan(theta_i) * np.tan(theta_r) * np.cos(phi_r) * script_n)

    p_d = rho_L / np.pi * N * W * np.cos(theta_r)

    p_s = (1 - Lambda) * 1 / (4 * np.cos(theta_i)) * F_ * P * G

    # integrated over the photodiode solid angle, the delta functions go away
    # in paper:
    # Lambda * F_ * G * 1 / (np.power(r, 2) * np.sin(theta_i)) * delta(theta_r - theta_i) * delta(phi_r)

    # assuming photodiode area is a slit with 10 times width as height
    photodiode_angle_radians = photodiode_angle * np.pi / 180.
    if np.abs(theta_r - theta_i) <= photodiode_angle_radians and np.abs(phi_r) <= 10 * photodiode_angle_radians:
        p_c = Lambda * F_ * G
        #print(theta_i * 180 / np.pi, theta_r * 180 / np.pi, photodiode_angle)
    else:
        p_c = 0

    return [p_d, p_s, p_c]


def BRIDF_diffuse(theta_r, phi_r, theta_i, n_0, polarization, photodiode_angle, parameters):
    return BRIDF_trio(theta_r, phi_r, theta_i, n_0, polarization, photodiode_angle, parameters)[0]


def BRIDF_specular_lobe(theta_r, phi_r, theta_i, n_0, polarization, photodiode_angle, parameters):
    return BRIDF_trio(theta_r, phi_r, theta_i, n_0, polarization, photodiode_angle, parameters)[1]


def BRIDF_specular_spike(theta_r, phi_r, theta_i, n_0, polarization, photodiode_angle, parameters):
    return BRIDF_trio(theta_r, phi_r, theta_i, n_0, polarization, photodiode_angle, parameters)[2]


def BRIDF(theta_r, phi_r, theta_i, n_0, polarization, photodiode_angle, parameters):
    trio = BRIDF_trio(theta_r, phi_r, theta_i, n_0, polarization, photodiode_angle, parameters)
    return trio[0] + trio[1] + trio[2]


def BRIDF_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees,
                  n_0, polarization, photodiode_angle, parameters):

    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        theta_r = np.pi * theta_r_in_degrees / 180
        return_array.append(BRIDF(theta_r, phi_r, theta_i, n_0, polarization, photodiode_angle, parameters))
    return return_array


def BRIDF_diffuse_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees,
                  n_0, polarization, photodiode_angle, parameters):
    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        theta_r = np.pi * theta_r_in_degrees / 180
        return_array.append(BRIDF_diffuse(theta_r, phi_r, theta_i, n_0, polarization, photodiode_angle, parameters))
    return return_array


def BRIDF_specular_lobe_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees,
                  n_0, polarization, photodiode_angle, parameters):
    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        theta_r = np.pi * theta_r_in_degrees / 180
        return_array.append(
            BRIDF_specular_lobe(theta_r, phi_r, theta_i, n_0, polarization, photodiode_angle, parameters))
    return return_array


def BRIDF_specular_spike_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees,
                  n_0, polarization, photodiode_angle, parameters):
    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        theta_r = np.pi * theta_r_in_degrees / 180
        return_array.append(
            BRIDF_specular_spike(theta_r, phi_r, theta_i, n_0, polarization, photodiode_angle, parameters))
    return return_array


# independent variables has the form
# [theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, photodiode_angle]
def unvectorized_fitter(independent_variables, log_rho_L, log_n_minus_one, log_K, log_gamma):
    theta_r = independent_variables[0] * np.pi / 180
    phi_r = independent_variables[1] * np.pi / 180
    theta_i = independent_variables[2] * np.pi / 180
    n_0 = independent_variables[3]
    polarization = independent_variables[4]
    photodiode_angle = independent_variables[5]
    rho_L = np.exp(log_rho_L)
    n = np.exp(log_n_minus_one) + 1
    K = np.exp(log_K)
    gamma = np.exp(log_gamma)
    parameters = [rho_L, n, K, gamma]
    return BRIDF(theta_r, phi_r, theta_i, n_0, polarization, photodiode_angle, parameters)


# takes arrays of independent variable lists
def fitter(independent_variables_array, log_rho_L, log_n_minus_one, log_K, log_gamma):
    arr = []
    for independent_variables in independent_variables_array:
        arr.append(unvectorized_fitter(independent_variables, log_rho_L, log_n_minus_one, log_K, log_gamma))
    return arr


# each point has form
# [theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity]
def fit_parameters(points):
    independent_variables_array = []
    intensity_array = []
    for point in points:
        independent_variables_array.append([point.theta_r_in_degrees, point.phi_r_in_degrees,
                                point.theta_i_in_degrees, point.n_0, point.polarization, point.photodiode_angular_width])
        intensity_array.append(point.intensity)
    # initial parameters are the ones found in the paper
    fit_params = scipy.optimize.curve_fit(fitter, independent_variables_array, intensity_array,
                                        p0=[np.log(0.5), np.log(1.5 - 1), np.log(0.4), np.log(0.05)])[0]
    return [np.exp(fit_params[0]), np.exp(fit_params)[1] + 1, np.exp(fit_params[2]), np.exp(fit_params[3])]
