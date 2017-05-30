import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

# convention
phi_r = 0


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
    # takes data = [[xdata1, ydata1, n_0_1, theta_i1], [...], ...]
    def __init__(self, data):
        self.data = data

    # BRIDF with n definitely greater than one and rho_L positive to be used for fitting
    # piecewise BRIDF of each test sequentially
    def fitter(self, theta_r_in_degrees_array, log_rho_L, log_n_minus_one, log_gamma):
        arr = []
        for data in self.data:
            arr += BRIDF(theta_r_in_degrees_array, self.data[3], self.data[2],
                        np.exp(log_rho_L), np.exp(log_n_minus_one) + 1, np.exp(log_gamma))
        return [x[0] + x[1] for x in arr]

# BRIDF as described in paper
# returns [[specular1, diffuse1], [specular2, diffuse2], ...]
def BRIDF(theta_r_in_degrees_array, theta_i_in_degrees, n_0, rho_L, n, gamma):
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    alpha = alpha_getter(gamma)
    for theta_r_in_degrees in theta_r_in_degrees_array:
        arr = BRIDF_helper(np.pi * theta_r_in_degrees / 180, theta_i, n_0, rho_L, n, gamma, alpha)
        return_array.append([arr[0], arr[1]])
    return return_array


# uses radians, has a random alpha already chosen
# returns [specular, diffuse]
def BRIDF_helper(theta_r, theta_i, n_0, rho_L, n, gamma, alpha):
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

    G = 1
    # temporary, G seems to always be zero
    """
    def G_prime(theta):
        return 2 / (1 + np.sqrt(1 + np.power(gamma * np.tan(theta), 2)))
    
    G = H(theta_i_prime - np.pi / 2) * H(theta_r_prime - np.pi / 2) * \
        G_prime(theta_i) * G_prime(theta_r)
    """

    F_ = F(theta_i, n / n_0)

    return [F_ * G * P_ / (4 * np.cos(theta_i)),
        rho_L / np.pi * W * (1 - A + B) * np.cos(theta_r)]


def fit_parameters(xdata, ydata, theta_i, n_0):

    fit = Fitter(theta_i, n_0)

    fit_params = scipy.optimize.curve_fit(fit.fitter, xdata, ydata,
                                          p0=[np.log(0.5), np.log(1.5 - 1), np.log(0.05)])[0]
    return [np.exp(fit_params[0]), np.exp(fit_params[1]) + 1, np.exp(fit_params[2])]


# used to fit parameters for multiple incident angles and n_0
# takes data = [[xdata1, ydata1, n_0_1, theta_i1], [...], ...]
# returns [rho_L, n, gamma]
# haven't tested this yet, waiting for more data
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


def plot_data(xdata, ydata, theta_i, n_0, title):
    fit_params = fit_parameters(xdata, ydata, theta_i, n_0)

    fitted_data = BRIDF(xdata, theta_i, n_0, fit_params[0], fit_params[1], fit_params[2])

    plt.plot(xdata, ydata, label="experimental")
    plt.plot(xdata, [x[0] for x in fitted_data], label="fitted specular")
    plt.plot(xdata, [x[1] for x in fitted_data], label="fitted diffuse")
    plt.plot(xdata, [x[0] + x[1] for x in fitted_data], label="fitted total")

    plt.title(title)
    plt.legend(bbox_to_anchor=(0.9, 0.9), bbox_transform=plt.gcf().transFigure)
    plt.ylabel("relative intensity")
    plt.xlabel("viewing angle (degrees)")
    string = "theta_i: " + str(theta_i) + "\nrho_L: " + str(fit_params[0]) + "\nn: " + \
        str(fit_params[1]) + "\ngamma: " + str(fit_params[2])
    plt.annotate(string, xy=(0.05, 0.8), xycoords='axes fraction')
    plt.show()

