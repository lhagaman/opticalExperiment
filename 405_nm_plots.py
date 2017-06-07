import TSTR_fit
import gaussian_fit
import semi_empirical_fit
import numpy as np
import matplotlib.pyplot as plt

read_data = True
plot_all_original_data = False
plot_all_original_data_with_paper_data_fits = False
plot_semi_empirical_with_all_points = True
plot_semi_empirical_with_all_points_adjusted = False



if read_data:

    dir = "405 nm Blue Laser diode reflectance measurements/"

    data_30_degree_1 = np.loadtxt(dir + "30 deg with blue laser diode run 2 after angle change 6-6.txt", skiprows=1)

    data_30_degree_2 = np.loadtxt(dir + "30 deg with blue laser diode run 2 after angle change.txt", skiprows=1)

    data_30_degree_3 = np.loadtxt(dir + "30 deg with blue laser diode.txt", skiprows=1)

    data_45_degree_1 = np.loadtxt(dir + "45 deg with blue laser diode 6-5 run 2 after changing angles.txt", skiprows=1)

    data_45_degree_2 = np.loadtxt(dir + "45 deg with blue laser diode 6-5.txt", skiprows=1)

    data_45_degree_3 = np.loadtxt(dir + "45 deg with blue laser diode.txt", skiprows=1)

    data_45_degree_4 = np.loadtxt(dir + "45 deg.txt", skiprows=1)

    data_30_degree_1_x = [np.rint(point[0]) for point in data_30_degree_1]
    data_30_degree_1_y = [point[1] for point in data_30_degree_1]

    data_30_degree_2_x = [np.rint(point[0]) for point in data_30_degree_2]
    data_30_degree_2_y = [point[1] for point in data_30_degree_2]

    data_30_degree_3_x = [np.rint(point[0]) for point in data_30_degree_3]
    data_30_degree_3_y = [point[1] for point in data_30_degree_3]

    data_45_degree_1_x = [np.rint(point[0]) for point in data_45_degree_1]
    data_45_degree_1_y = [point[1] for point in data_45_degree_1]

    data_45_degree_2_x = [np.rint(point[0]) for point in data_45_degree_2]
    data_45_degree_2_y = [point[1] for point in data_45_degree_2]

    data_45_degree_3_x = [np.rint(point[0]) for point in data_45_degree_3]
    data_45_degree_3_y = [point[1] for point in data_45_degree_3]

    data_45_degree_4_x = [np.rint(point[0]) for point in data_45_degree_4]
    data_45_degree_4_y = [point[1] for point in data_45_degree_4]

if plot_all_original_data:

    f, axarr = plt.subplots(2, sharex=True)

    axarr[0].plot(data_30_degree_1_x, data_30_degree_1_y, label="run 2 after angle change 6-6")
    axarr[0].plot(data_30_degree_2_x, data_30_degree_2_y, label="run 2 after angle change")
    axarr[0].plot(data_30_degree_3_x, data_30_degree_3_y, label="before angle change")

    axarr[0].set_title("30 deg with blue laser diode")
    axarr[0].legend()

    axarr[1].plot(data_45_degree_1_x, data_45_degree_1_y, label="6-5 run 2 after changing angles")
    axarr[1].plot(data_45_degree_2_x, data_45_degree_2_y, label="6-5")
    axarr[1].plot(data_45_degree_3_x, data_45_degree_3_y, label="run 2")
    axarr[1].plot(data_45_degree_4_x, data_45_degree_4_y, label="run 1")

    axarr[1].set_title("45 deg with blue laser diode")
    axarr[1].legend()

    plt.ylabel("relative intensity")
    plt.xlabel("viewing angle (degrees)")
    plt.show()

if plot_all_original_data_with_paper_data_fits:

    fig, axarr = plt.subplots(2, 2)
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.ylabel("relative intensity")
    plt.xlabel("viewing angle (degrees)")
    ax = plt.gca()
    ax.yaxis.set_label_coords(-0.1, 0.5)

    axarr[0, 0].plot(data_30_degree_1_x, data_30_degree_1_y, label="run 2 after angle change 6-6")
    axarr[0, 0].plot(data_30_degree_2_x, data_30_degree_2_y, label="run 2 after angle change")
    axarr[0, 0].plot(data_30_degree_3_x, data_30_degree_3_y, label="before angle change")
    axarr[0, 0].set_title("30 deg with blue laser diode")
    axarr[0, 0].legend(prop={'size':6})

    axarr[1, 0].plot(data_45_degree_1_x, data_45_degree_1_y, label="6-5 run 2 after changing angles")
    axarr[1, 0].plot(data_45_degree_2_x, data_45_degree_2_y, label="6-5")
    axarr[1, 0].plot(data_45_degree_3_x, data_45_degree_3_y, label="run 2")
    axarr[1, 0].plot(data_45_degree_4_x, data_45_degree_4_y, label="run 1")
    axarr[1, 0].set_title("45 deg with blue laser diode")
    axarr[1, 0].legend(prop={'size':6})

    parameters = [0.7, 1.5, 3, 0.055]
    axarr[0, 1].plot(data_30_degree_1_x, semi_empirical_fit.BRIDF_plotter(
        data_30_degree_1_x, 0, 30, 1, 0, parameters),
             label="Semi-Empirical\nusing parameters from paper")
    parameters = [0.1768, 40.59, 48.99]
    y_fit = gaussian_fit.BRIDF_plotter(data_30_degree_1_x, 30, parameters)
    y_fit = [0.015 * y for y in y_fit]
    axarr[0, 1].plot(data_30_degree_1_x, y_fit, label="Gaussian\nusing parameters from paper \n(hand-chosen normalization)")
    axarr[0, 1].legend(prop={'size':6})

    parameters = [0.7, 1.5, 3, 0.055]
    axarr[1, 1].plot(data_45_degree_1_x, semi_empirical_fit.BRIDF_plotter(
        data_45_degree_1_x, 0, 45, 1, 0, parameters),
                     label="Semi-Empirical\nusing parameters from paper")
    parameters = [0.1469, 90, 50.18]
    y_fit = gaussian_fit.BRIDF_plotter(data_45_degree_1_x, 45, parameters)
    y_fit = [0.015 * y for y in y_fit]
    axarr[1, 1].plot(data_45_degree_1_x, y_fit,
                     label="Gaussian\nusing parameters from paper \n(hand-chosen normalization)")
    axarr[1, 1].legend(prop={'size':6})

    plt.show()


# each point has form
# [theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity]
def make_points(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity_array):
    points = []
    numpoints = len(theta_r_in_degrees_array)
    for i in range(numpoints):
        theta_r_in_degrees = theta_r_in_degrees_array[i]
        # try to somewhat account fot the fact that we're not using the sr^-1 units we should be using
        intensity = 100 * intensity_array[i]
        points.append([theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees,
                       n_0, polarization, intensity])
    return points

points_30_degree_1 = make_points(data_30_degree_1_x, 0, 30, 1, 0, data_30_degree_1_y)
points_30_degree_2 = make_points(data_30_degree_2_x, 0, 30, 1, 0, data_30_degree_2_y)
points_30_degree_3 = make_points(data_30_degree_3_x, 0, 30, 1, 0, data_30_degree_3_y)
points_45_degree_1 = make_points(data_45_degree_1_x, 0, 45, 1, 0, data_45_degree_1_y)
points_45_degree_2 = make_points(data_45_degree_2_x, 0, 45, 1, 0, data_45_degree_2_y)
points_45_degree_3 = make_points(data_45_degree_3_x, 0, 45, 1, 0, data_45_degree_3_y)
points_45_degree_4 = make_points(data_45_degree_4_x, 0, 45, 1, 0, data_45_degree_4_y)

all_points = points_30_degree_1 + points_30_degree_2 + points_30_degree_3 + points_45_degree_1 + \
             points_45_degree_2 + points_45_degree_3 + points_45_degree_4

# here I manually adjust some input angles to see how it affects the fit
# because I suspect some have large errors
adjusted_points_30_degree_1 = make_points(data_30_degree_1_x, 0, 32, 1, 0, data_30_degree_1_y)
adjusted_points_30_degree_2 = make_points(data_30_degree_2_x, 0, 32, 1, 0, data_30_degree_2_y)
adjusted_points_30_degree_3 = make_points(data_30_degree_3_x, 0, 41, 1, 0, data_30_degree_3_y)
adjusted_points_45_degree_2 = make_points(data_45_degree_2_x, 0, 50, 1, 0, data_45_degree_2_y)
adjusted_points_45_degree_3 = make_points(data_45_degree_3_x, 0, 51, 1, 0, data_45_degree_3_y)
adjusted_points_45_degree_4 = make_points(data_45_degree_4_x, 0, 48, 1, 0, data_45_degree_4_y)

all_points_adjusted = adjusted_points_30_degree_1 + adjusted_points_30_degree_2 + adjusted_points_30_degree_3 + \
                      points_45_degree_1 + adjusted_points_45_degree_2 + adjusted_points_45_degree_3 + \
                      adjusted_points_45_degree_4

def plot_experimental_data(points, title):
    x_data = [point[0] for point in points]
    ydata = [point[5] for point in points]
    plt.plot(x_data, ydata)
    plt.title(title)
    plt.show()


# implicitly assumes only theta_r, theta_i, and intensity vary, otherwise will get weird results
# only does gaussian fit if there's only one theta_i
# (because the fit parameters are functions of theta_i, more complicated)
def plot_with_semi_empirical_and_gaussian_fits(points, title):
    one_pass_x_data = []
    theta_i_list = []
    points_by_theta_i = []
    for point in points:
        if point[0] not in one_pass_x_data:
            one_pass_x_data.append(point[0])
        if point[2] not in theta_i_list:
            theta_i_list.append(point[2])
            points_by_theta_i.append([point])
        else:
            index = 0.5 # to error if not overridden
            for i in range(len(theta_i_list)):
                theta_i = theta_i_list[i]
                if np.round(point[2], 2) == theta_i:
                    index = i
            points_by_theta_i[index].append(point)

    plot_gaussian = len(theta_i_list) == 1

    phi_r_in_degrees = points[0][1]
    n_0 = points[0][3]
    polarization = points[0][4]

    semi_empirical_parameters = semi_empirical_fit.fit_parameters(points)
    if plot_gaussian:
        gaussian_parameters = gaussian_fit.fit_parameters(points)

    for i in range(len(theta_i_list)):

        theta_i_in_degrees = theta_i_list[i]
        points = points_by_theta_i[i]

        x_data = [point[0] for point in points]
        y_data = [point[5] for point in points]

        semi_empirical_y = semi_empirical_fit.BRIDF_plotter(one_pass_x_data,
                                phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, semi_empirical_parameters)
        if plot_gaussian:
            gaussian_y = gaussian_fit.BRIDF_plotter(one_pass_x_data, theta_i_in_degrees, gaussian_parameters)

        plt.figure()
        plt.title(title)
        plt.scatter(x_data, y_data, marker="x", color="g", label="experimental")
        plt.plot(one_pass_x_data, semi_empirical_y, label="Semi-Empirical Fit")
        if plot_gaussian:
            plt.plot(one_pass_x_data, gaussian_y, label="Gaussian Fit")

        rho_L = semi_empirical_parameters[0]
        n = semi_empirical_parameters[1]
        K = semi_empirical_parameters[2]
        gamma = semi_empirical_parameters[3]

        if plot_gaussian:
            sigma = gaussian_parameters[0]
            R_1 = gaussian_parameters[1]
            R_2 = gaussian_parameters[2]

        string = "theta_i: " + str(theta_i_in_degrees) + "\n\nrho_L: " + str(rho_L) + "\nn: " + \
                 str(n) + "\nK: " + str(K) + "\ngamma: " + str(gamma)
        if plot_gaussian:
            string += "\n\nsigma: " + str(sigma) + "\nR_1: " + str(R_1) + "\nR_2: " + str(R_2)
        plt.legend()
        plt.xlabel("viewing angle (degrees)")
        plt.ylabel("relative intensity")
        plt.annotate(string, xy=(0.05, 0.7), xycoords='axes fraction', size=6)

    plt.show()

# these fit poorly
# plot_with_semi_empirical_and_gaussian_fits(points_30_degree_1,
#                                           "30 deg with blue laser diode run 2 after angle change 6-6.txt")
# plot_with_semi_empirical_and_gaussian_fits(points_45_degree_2, "45 deg with blue laser diode 6-5")
# plot_with_semi_empirical_and_gaussian_fits(points_45_degree_3, "45 deg with blue laser diode")
# plot_with_semi_empirical_and_gaussian_fits(points_45_degree_4, "45 deg")


# these fit well
# plot_with_semi_empirical_and_gaussian_fits(points_45_degree_1,
#                                           "45 deg with blue laser diode 6-5 run 2 after changing angles")
# plot_with_semi_empirical_and_gaussian_fits(adjusted_points_45_degree_2,
#                                           "45 deg with blue laser diode 6-5, adjusted theta_i")
# plot_with_semi_empirical_and_gaussian_fits(adjusted_points_45_degree_3,
#                                           "45 deg with blue laser diode, adjusted theta_i")
# plot_with_semi_empirical_and_gaussian_fits(adjusted_points_45_degree_4, "45 deg, adjusted theta_i")

# plot_with_semi_empirical_and_gaussian_fits(all_points, "All Points")
# plot_with_semi_empirical_and_gaussian_fits(all_points_adjusted, "All Points With Manually Adjusted theta_i Values")
