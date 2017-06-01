import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt


def F(theta, n_frac):
    theta_t = np.arcsin(1 / n_frac * np.sin(theta))
    return 0.5 * np.power((np.sin(theta - theta_t) / np.sin(theta + theta_t)), 2) * \
        (1 + np.power((np.cos(theta + theta_t) / np.cos(theta - theta_t)), 2))


# heaviside step function
def H(x):
    return 0.5 * (np.sign(x) + 1)


# BRIDF as described in paper
# takes an array, easy to plot
# returns [[specular1, diffuse1], [specular2, diffuse2], ...]
def BRIDF_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, rho_L, n, gamma):
    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        arr = BRIDF(np.pi * theta_r_in_degrees / 180, phi_r, theta_i, n_0, rho_L, n, gamma)
        return_array.append([arr[0], arr[1]])
    return return_array


# uses radians, takes a single point
# returns [specular, diffuse]
def BRIDF(theta_r, phi_r, theta_i, n_0, rho_L, n, gamma):
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

    # There was a typo here, I changed by moving parentheses and taking a reciprocal
    W = (1 - F(theta_i, n / n_0)) * (1 - F(np.arcsin(n_0 / n * np.sin(theta_r)), n / n_0))

    A = 0.5 * np.power(gamma, 2) / (np.power(gamma, 2) + 0.92)

    theta_m = min(theta_i, theta_r)

    theta_M = max(theta_i, theta_r)

    B = 0.5 * np.power(gamma, 2) / (np.power(gamma, 2) + 0.25) * \
        H(np.cos(phi_r)) * np.cos(phi_r) * np.sin(theta_M) * np.tan(theta_m)

    def G_prime(theta):
        return 2 / (1 + np.sqrt(1 + np.power(gamma * np.tan(theta), 2)))

    # this has negative of the inside of the H functions from the paper, I think it was a typo
    G = H(np.pi / 2 - theta_i_prime) * H(np.pi / 2 - theta_r_prime) * \
        G_prime(theta_i) * G_prime(theta_r)

    F_ = F(theta_i, n / n_0)

    return [F_ * G * P_ / (4 * np.cos(theta_i)),
        rho_L / np.pi * W * (1 - A + B) * np.cos(theta_r)]


# independent variables has the form [theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0]
def unvectorized_fitter(independent_variables, log_rho_L, log_n_minus_one, log_gamma):
    theta_r = independent_variables[0] * np.pi / 180
    phi_r = independent_variables[1] * np.pi / 180
    theta_i = independent_variables[2] * np.pi / 180
    n_0 = independent_variables[3]
    rho_L = np.exp(log_rho_L)
    n = np.exp(log_n_minus_one) + 1
    gamma = np.exp(log_gamma)
    pair = BRIDF(theta_r, phi_r, theta_i, n_0, rho_L, n, gamma)
    return pair[0] + pair[1]


# takes arrays of independent variable lists
def fitter(independent_variables_array, log_rho_L, log_n_minus_one, log_gamma):
    arr = []
    for independent_variables in independent_variables_array:
        arr.append(unvectorized_fitter(independent_variables, log_rho_L, log_n_minus_one, log_gamma))
    return arr


# each point has form [theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, intensity]
def fit_parameters(points):
    independent_variables_array = []
    intensity_array = []
    for point in points:
        independent_variables_array.append(point[:4])
        intensity_array.append(point[4])
    # initial parameters are the ones found in the paper
    fit_params = scipy.optimize.curve_fit(fitter, independent_variables_array, intensity_array,
                                          p0=[np.log(0.5), np.log(1.5 - 1), np.log(0.05)])[0]
    return [np.exp(fit_params[0]), np.exp(fit_params)[1] + 1, np.exp(fit_params[2])]


# works for plotting curves with constant phi_r, theta_i, and n_0
def plot_data_with_fit(points, title):

    phi_r_in_degrees = points[0][1]
    theta_i_in_degrees = points[0][2]
    n_0 = points[0][3]

    parameters = fit_parameters(points)

    rho_L = parameters[0]
    n = parameters[1]
    gamma = parameters[2]

    theta_r_in_degrees_array = []
    intensity_array = []
    for point in points:
        theta_r_in_degrees_array.append(point[0])
        intensity_array.append(point[4])

    fitted_data = BRIDF_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, rho_L, n, gamma)

    plt.plot(theta_r_in_degrees_array, intensity_array, label="experimental")
    plt.plot(theta_r_in_degrees_array, [x[0] for x in fitted_data], label="fitted specular")
    plt.plot(theta_r_in_degrees_array, [x[1] for x in fitted_data], label="fitted diffuse")
    plt.plot(theta_r_in_degrees_array, [x[0] + x[1] for x in fitted_data], label="fitted total")

    plt.title(title)
    plt.legend(bbox_to_anchor=(0.9, 0.9), bbox_transform=plt.gcf().transFigure)
    plt.ylabel("relative intensity")
    plt.xlabel("viewing angle (degrees)")
    string = "theta_i: " + str(theta_i_in_degrees) + "\nrho_L: " + str(rho_L) + "\nn: " + \
             str(n) + "\ngamma: " + str(gamma)
    plt.annotate(string, xy=(0.05, 0.8), xycoords='axes fraction')
    plt.show()
