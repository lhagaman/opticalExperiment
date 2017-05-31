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


class RandomAlphaChooser:

    def __init__(self, gamma):

        def P(alpha_):
            return np.power(gamma, 2) / \
                (np.pi * np.power(np.cos(alpha_), 4) *
                    np.power(np.power(gamma, 2) + np.power(np.tan(alpha_), 2), 2))

        def alpha_pdf(alpha_):
            return P(alpha_) * np.sin(alpha_)

        random_alpha_values = []

        even_alpha_values = np.linspace(0, np.pi/2, 100)

        if gamma < 0.01:
            self.alphas = [0]
        else:

            for alpha in even_alpha_values:
                for i in range(np.rint(100 * alpha_pdf(alpha)).astype(int)):
                    random_alpha_values.append(alpha)

            self.alphas = random_alpha_values

    def get_alpha(self):
        return np.random.choice(self.alphas)


def alpha_getter(gamma):

    chooser = RandomAlphaChooser(gamma)

    N = 1000
    sum = 0
    for i in range(N):
        sum += chooser.get_alpha()

    return sum / N

class Fitter:
    def __init__(self, theta_i, n_0):
        self.theta_i = theta_i
        self.n_0 = n_0

    # BRIDF with n definitely greater than one and rho_L positive to be used for fitting
    def fitter(self, theta_r_in_degrees_array, log_rho_L, log_n_minus_one, log_gamma):
        arr = BRIDF(theta_r_in_degrees_array, self.theta_i, self.n_0,
                    np.exp(log_rho_L), np.exp(log_n_minus_one) + 1, np.exp(log_gamma))
        return [x[0] + x[1] for x in arr]


class MultiFitter:
    # takes data = [[xdata1, ydata1, n_0_1, theta_i1, phi_r1], [...], ...]
    def __init__(self, data):
        self.data = data

    # BRIDF with n definitely greater than one and rho_L positive to be used for fitting
    # piecewise BRIDF of each test sequentially
    def fitter(self, theta_r_in_degrees_array, log_rho_L, log_n_minus_one, log_gamma):
        arr = []
        for data in self.data:
            arr += BRIDF(theta_r_in_degrees_array, data[4], data[3], data[2],
                        np.exp(log_rho_L), np.exp(log_n_minus_one) + 1, np.exp(log_gamma))
        return [x[0] + x[1] for x in arr]


# BRIDF as described in paper
# returns [[specular1, diffuse1], [specular2, diffuse2], ...]
def BRIDF(theta_r_in_degrees_array, phi_r, theta_i_in_degrees, n_0, rho_L, n, gamma):
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        arr = BRIDF_helper(np.pi * theta_r_in_degrees / 180, phi_r, theta_i, n_0, rho_L, n, gamma)
        return_array.append([arr[0], arr[1]])
    return return_array


# uses radians, has a random alpha already chosen
# returns [specular, diffuse]
def BRIDF_helper(theta_r, phi_r, theta_i, n_0, rho_L, n, gamma):
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


# used to fit parameters for multiple incident angles and n_0
# takes data = [[xdata1, ydata1, n_0_1, theta_i1, phi_r1], [...], ...]
# returns [rho_L, n, gamma]
# haven't tested this with multiple datasets yet, waiting for more data
def multi_fit_parameters(data):
    fit = MultiFitter(data)

    totalxdata = []
    totalydata = []

    for x in data:
        totalxdata += x[0]
        totalydata += x[1]

    fit_params = scipy.optimize.curve_fit(fit.fitter, totalxdata, totalydata,
                                          p0=[np.log(0.5), np.log(1.5 - 1), np.log(0.05)])[0]
    return [np.exp(fit_params[0]), np.exp(fit_params[1]) + 1, np.exp(fit_params[2])]


# for plotting just one theta_i and n_0
def plot_data(xdata, ydata, n_0, theta_i, phi_r, title):

    data = [[xdata, ydata, n_0, theta_i, phi_r]]

    totalxdata = []
    totalydata = []

    for x in data:
        totalxdata += x[0]
        totalydata += x[1]

    fit_params = multi_fit_parameters(data)

    fitted_data = BRIDF(xdata, phi_r, theta_i, n_0, fit_params[0], fit_params[1], fit_params[2])

    plt.plot(totalxdata, totalydata, label="experimental")
    plt.plot(totalxdata, [x[0] for x in fitted_data], label="fitted specular")
    plt.plot(totalxdata, [x[1] for x in fitted_data], label="fitted diffuse")
    plt.plot(totalxdata, [x[0] + x[1] for x in fitted_data], label="fitted total")

    plt.title(title)
    plt.legend(bbox_to_anchor=(0.9, 0.9), bbox_transform=plt.gcf().transFigure)
    plt.ylabel("relative intensity")
    plt.xlabel("viewing angle (degrees)")
    string = "theta_i: " + str(theta_i) + "\nrho_L: " + str(fit_params[0]) + "\nn: " + \
        str(fit_params[1]) + "\ngamma: " + str(fit_params[2])
    plt.annotate(string, xy=(0.05, 0.8), xycoords='axes fraction')
    plt.show()

