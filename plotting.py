import semi_empirical_fit, TSTR_fit, gaussian_fit, large_gas_layer_fit, partial_gas_layer_fit, reff_polynomial_fit
import numpy as np
import matplotlib.pyplot as plt
from Point import Point


# implicitly assumes only theta_r, theta_i, and intensity vary, otherwise will get weird results
# only does gaussian fit if there's only one theta_i
# (because the fit parameters are functions of theta_i, more complicated)
def plot_with_semi_empirical_and_gaussian_fits(points):
    one_pass_x_data = []
    run_name_list = []
    points_by_run_name = []
    for point in points:
        if point.theta_r_in_degrees not in one_pass_x_data:
            one_pass_x_data.append(point.theta_r_in_degrees)
        if point.run_name not in run_name_list:
            run_name_list.append(point.run_name)
            points_by_run_name.append([point])
        else:
            index = 0.5  # to error if not overridden
            for i in range(len(run_name_list)):
                run_name = run_name_list[i]
                if point.run_name == run_name:
                    index = i
            points_by_run_name[index].append(point)

    plot_gaussian = len(run_name_list) == 1

    phi_r_in_degrees = points[0].phi_r_in_degrees
    n_0 = points[0].n_0
    polarization = points[0].polarization
    photodiode_solid_angle = points[0].photodiode_solid_angle

    semi_empirical_parameters = semi_empirical_fit.fit_parameters(points)
    TSTR_parameters = TSTR_fit.fit_parameters(points)
    if plot_gaussian:
        gaussian_parameters = gaussian_fit.fit_parameters(points)

    for i in range(len(run_name_list)):

        run_name = run_name_list[i]
        points = points_by_run_name[i]
        theta_i_in_degrees = points[0].theta_i_in_degrees

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]
        max_y = max(y_data)

        semi_empirical_y = semi_empirical_fit.BRIDF_plotter(one_pass_x_data,
                                                            phi_r_in_degrees, theta_i_in_degrees, n_0, polarization,
                                                            photodiode_solid_angle,
                                                            semi_empirical_parameters)
        TSTR_y = TSTR_fit.BRIDF_plotter(one_pass_x_data,
                                        phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, TSTR_parameters)
        if plot_gaussian:
            gaussian_y = gaussian_fit.BRIDF_plotter(one_pass_x_data, theta_i_in_degrees, gaussian_parameters)

        plt.figure()
        plt.title(run_name)
        plt.scatter(x_data, y_data, marker="x", color="r", label="experimental")
        plt.plot(one_pass_x_data, semi_empirical_y, label="Semi-Empirical Fit")
        plt.plot(one_pass_x_data, TSTR_y, label="TSTR Fit")
        if plot_gaussian:
            plt.plot(one_pass_x_data, gaussian_y, label="Gaussian Fit")

        rho_L_ = semi_empirical_parameters[0]
        n_ = semi_empirical_parameters[1]
        K_ = semi_empirical_parameters[2]
        gamma_ = semi_empirical_parameters[3]
        rho_L = TSTR_parameters[0]
        n = TSTR_parameters[1]
        gamma = TSTR_parameters[2]

        if plot_gaussian:
            sigma = gaussian_parameters[0]
            R_1 = gaussian_parameters[1]
            R_2 = gaussian_parameters[2]

        string = "theta_i: " + str(theta_i_in_degrees) + "\n\nSemi-Empirical Parameters:\nrho_L: " + str(rho_L_) + \
                 "\nn: " + str(n_) + "\nK: " + str(K_) + "\ngamma: " + str(gamma_) + "\n\nTSTR Parameters:\nrho_L: " + \
                 str(rho_L) + "\nn: " + str(n) + "\ngamma: " + str(gamma)
        if plot_gaussian:
            string += "\n\nsigma: " + str(sigma) + "\nR_1: " + str(R_1) + "\nR_2: " + str(R_2)
        axes = plt.gca()
        axes.set_ylim([0, 1.2 * max_y])
        plt.legend()
        plt.xlabel("viewing angle (degrees)")
        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.annotate(string, xy=(0.05, 0.6), xycoords='axes fraction', size=6)

    plt.show()


def plot_experimental_data(points, title):
    x_data = [point[0] for point in points]
    ydata = [point[5] for point in points]
    plt.plot(x_data, ydata)
    plt.title(title)
    plt.show()


def plot_semi_empirical_components(points):
    one_pass_x_data = []
    run_name_list = []
    points_by_run_name = []
    for point in points:
        if point.theta_r_in_degrees not in one_pass_x_data:
            one_pass_x_data.append(point.theta_r_in_degrees)
        if point.run_name not in run_name_list:
            run_name_list.append(point.run_name)
            points_by_run_name.append([point])
        else:
            index = 0.5  # to error if not overridden
            for i in range(len(run_name_list)):
                run_name = run_name_list[i]
                if point.run_name == run_name:
                    index = i
            points_by_run_name[index].append(point)

    phi_r_in_degrees = points[0].phi_r_in_degrees
    n_0 = points[0].n_0
    polarization = points[0].polarization
    photodiode_solid_angle = points[0].photodiode_solid_angle

    semi_empirical_parameters = semi_empirical_fit.fit_parameters(points)

    for i in range(len(run_name_list)):

        run_name = run_name_list[i]
        points = points_by_run_name[i]
        theta_i_in_degrees = points[0].theta_i_in_degrees

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]
        max_y = max(y_data)

        semi_empirical_y_diffuse = semi_empirical_fit.BRIDF_diffuse_plotter(one_pass_x_data, phi_r_in_degrees,
                                                                            theta_i_in_degrees, n_0, polarization,
                                                                            photodiode_solid_angle,
                                                                            semi_empirical_parameters)
        semi_empirical_y_specular_lobe = semi_empirical_fit.BRIDF_specular_lobe_plotter(one_pass_x_data,
                                                                                        phi_r_in_degrees,
                                                                                        theta_i_in_degrees, n_0,
                                                                                        polarization,
                                                                                        photodiode_solid_angle,
                                                                                        semi_empirical_parameters)
        semi_empirical_y_specular_spike = semi_empirical_fit.BRIDF_specular_spike_plotter(one_pass_x_data,
                                                                                          phi_r_in_degrees,
                                                                                          theta_i_in_degrees, n_0,
                                                                                          polarization,
                                                                                          photodiode_solid_angle,
                                                                                          semi_empirical_parameters)
        semi_empirical_y_total = semi_empirical_fit.BRIDF_plotter(one_pass_x_data, phi_r_in_degrees, theta_i_in_degrees,
                                                                  n_0, polarization, photodiode_solid_angle,
                                                                  semi_empirical_parameters)

        plt.figure()
        plt.title(run_name + "\nSemi-Empirical Fit Components")
        plt.scatter(x_data, y_data, marker="x", color="r", label="experimental")
        plt.plot(one_pass_x_data, semi_empirical_y_diffuse, label="diffuse")
        plt.plot(one_pass_x_data, semi_empirical_y_specular_lobe, label="specular_lobe")
        plt.plot(one_pass_x_data, semi_empirical_y_specular_spike, label="specular_spike")
        plt.plot(one_pass_x_data, semi_empirical_y_total, label="total")

        rho_L_ = semi_empirical_parameters[0]
        n_ = semi_empirical_parameters[1]
        K_ = semi_empirical_parameters[2]
        gamma_ = semi_empirical_parameters[3]

        string = "theta_i: " + str(theta_i_in_degrees) + "\n\nSemi-Empirical Parameters:\nrho_L: " + str(rho_L_) + \
                 "\nn: " + str(n_) + "\nK: " + str(K_) + "\ngamma: " + str(gamma_)

        axes = plt.gca()
        axes.set_ylim([0, 1.2 * max_y])
        plt.legend()
        plt.xlabel("viewing angle (degrees)")
        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.annotate(string, xy=(0.05, 0.8), xycoords='axes fraction', size=6)

    plt.show()


def plot_large_gas_layer_semi_empirical(points):
    one_pass_x_data = []
    run_name_list = []
    points_by_run_name = []
    for point in points:
        if point.theta_r_in_degrees not in one_pass_x_data:
            one_pass_x_data.append(point.theta_r_in_degrees)
        if point.run_name not in run_name_list:
            run_name_list.append(point.run_name)
            points_by_run_name.append([point])
        else:
            index = 0.5  # to error if not overridden
            for i in range(len(run_name_list)):
                run_name = run_name_list[i]
                if point.run_name == run_name:
                    index = i
            points_by_run_name[index].append(point)

    phi_r_in_degrees = points[0].phi_r_in_degrees
    n_gas = points[0].n_0
    n_liquid = 1.69
    polarization = points[0].polarization
    photodiode_solid_angle = points[0].photodiode_solid_angle

    semi_empirical_parameters = large_gas_layer_fit.fit_parameters_semi_empirical(points)

    for i in range(len(run_name_list)):
        run_name = run_name_list[i]
        points = points_by_run_name[i]
        theta_i_in_degrees = points[0].theta_i_in_degrees

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]
        max_y = max(y_data)

        fit_y = large_gas_layer_fit.BRIDF_plotter(semi_empirical_fit.BRIDF, one_pass_x_data, phi_r_in_degrees, theta_i_in_degrees,
                                                                  n_gas, n_liquid, polarization, photodiode_solid_angle,
                                                                  semi_empirical_parameters)

        plt.figure()
        plt.title(run_name + "\nSemi-Empirical Fit Components")
        plt.scatter(x_data, y_data, marker="x", color="r", label="experimental")
        plt.plot(one_pass_x_data, fit_y, label="fit")

        rho_L_ = semi_empirical_parameters[0]
        n_ = semi_empirical_parameters[1]
        K_ = semi_empirical_parameters[2]
        gamma_ = semi_empirical_parameters[3]

        string = "theta_i: " + str(theta_i_in_degrees) + "\nn_liquid: 1.69 (chosen)" + \
                 "\n\nSemi-Empirical Parameters:\nrho_L: " + str(rho_L_) + \
                 "\nn_gas (fitted): " + str(n_) + "\nK: " + str(K_) + "\ngamma: " + str(gamma_)

        axes = plt.gca()
        axes.set_ylim([0, 1.2 * max_y])
        plt.legend()
        plt.xlabel("viewing angle (degrees)")
        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.annotate(string, xy=(0.05, 0.7), xycoords='axes fraction', size=6)

    plt.show()


def plot_large_gas_layer_gaussian(points):
    one_pass_x_data = []
    run_name_list = []
    points_by_run_name = []
    for point in points:
        if point.theta_r_in_degrees not in one_pass_x_data:
            one_pass_x_data.append(point.theta_r_in_degrees)
        if point.run_name not in run_name_list:
            run_name_list.append(point.run_name)
            points_by_run_name.append([point])
        else:
            index = 0.5  # to error if not overridden
            for i in range(len(run_name_list)):
                run_name = run_name_list[i]
                if point.run_name == run_name:
                    index = i
            points_by_run_name[index].append(point)

    phi_r_in_degrees = points[0].phi_r_in_degrees
    n_gas = 1
    n_liquid = 1.69
    polarization = points[0].polarization
    photodiode_solid_angle = points[0].photodiode_solid_angle

    gaussian_parameters = large_gas_layer_fit.fit_parameters_gaussian(points)

    for i in range(len(run_name_list)):
        run_name = run_name_list[i]
        points = points_by_run_name[i]
        theta_i_in_degrees = points[0].theta_i_in_degrees

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]
        max_y = max(y_data)

        fit_y = large_gas_layer_fit.BRIDF_plotter(gaussian_fit.BRIDF_all_parameters, one_pass_x_data, phi_r_in_degrees,
                                                  theta_i_in_degrees, n_gas, n_liquid, polarization,
                                                  photodiode_solid_angle, gaussian_parameters)

        plt.figure()
        plt.title(run_name + "\nLarge Gas Layer / Gaussian Fit")
        plt.scatter(x_data, y_data, marker="x", color="r", label="experimental")
        plt.plot(one_pass_x_data, fit_y, label="fit")

        sigma = gaussian_parameters[0]
        R_1 = gaussian_parameters[1]
        R_2 = gaussian_parameters[2]

        string = "theta_i: " + str(theta_i_in_degrees) + "\nn_liquid: 1.69 (chosen)" + "\nn_gas: 1 (chosen)" + \
                 "\n\nGaussian Parameters:\nsigma: " + str(sigma) + \
                 "\nR_1: " + str(R_1) + "\nR_2: " + str(R_2)

        axes = plt.gca()
        axes.set_ylim([0, 1.2 * max_y])
        plt.legend()
        plt.xlabel("viewing angle (degrees)")
        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.annotate(string, xy=(0.05, 0.7), xycoords='axes fraction', size=6)

    plt.show()


def plot_partial_gas_layer_gaussian(points_with_liquid, points_without_liquid, title):
    one_pass_x_data = []
    run_name_list = []
    points_by_run_name = []
    points = points_with_liquid
    for point in points:
        if point.theta_r_in_degrees not in one_pass_x_data:
            one_pass_x_data.append(point.theta_r_in_degrees)
        if point.run_name not in run_name_list:
            run_name_list.append(point.run_name)
            points_by_run_name.append([point])
        else:
            index = 0.5  # to error if not overridden
            for i in range(len(run_name_list)):
                run_name = run_name_list[i]
                if point.run_name == run_name:
                    index = i
            points_by_run_name[index].append(point)

    phi_r_in_degrees = points[0].phi_r_in_degrees
    n_gas = 1
    n_liquid = 1.69
    polarization = points[0].polarization
    photodiode_solid_angle = points[0].photodiode_solid_angle

    parameters = partial_gas_layer_fit.fit_parameters_gaussian(points)

    sigma = parameters[0]
    R_1 = parameters[1]
    R_2 = parameters[2]
    x = parameters[3]

    for i in range(len(run_name_list)):
        run_name = run_name_list[i]
        points = points_by_run_name[i]
        theta_i_in_degrees = points[0].theta_i_in_degrees

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]

        fit_y = partial_gas_layer_fit.BRIDF_plotter(gaussian_fit.BRIDF_all_parameters, one_pass_x_data, phi_r_in_degrees,
                                                  theta_i_in_degrees, n_gas, n_liquid, polarization,
                                                  photodiode_solid_angle, [sigma, R_1, R_2], x)

        # this is the fit if there were no liquid
        fit_y_gas = gaussian_fit.BRIDF_plotter(one_pass_x_data, theta_i_in_degrees, [sigma, R_1, R_2])

        plt.figure()
        plt.title(run_name + "\nPartial Gas Layer / Gaussian Fit" + title)
        plt.scatter(x_data, y_data, marker="x", color="r", label="experimental with liquid")
        if points_without_liquid:
            x_data_ = [point.theta_r_in_degrees for point in points_without_liquid]
            y_data_ = [point.intensity for point in points_without_liquid]
            plt.scatter(x_data_, y_data_, marker="x", color="b", label="experimental without liquid")
        plt.plot(one_pass_x_data, fit_y, label="partial gas layer fit", color="r")
        plt.plot(one_pass_x_data, fit_y_gas,
                 label="gaussian fit if there\nwere no liquid (predicted from data with liquid)", color="b")

        string = "theta_i: " + str(theta_i_in_degrees) + "\nn_liquid: " + str(n_liquid) + "\nn_gas: " + str(n_gas) + \
                 "\n\nx: " + str(x) + "\nsigma: " + str(sigma) + \
                 "\nR_1: " + str(R_1) + "\nR_2: " + str(R_2)

        max_y = max(max(y_data), max(y_data_), max(fit_y), max(fit_y_gas))
        axes = plt.gca()
        axes.set_ylim([0, 1.2 * max_y])
        plt.legend()
        plt.xlabel("viewing angle (degrees)")
        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.annotate(string, xy=(0.05, 0.7), xycoords='axes fraction', size=6)


def plot_with_semi_empirical_TSTR_gaussian_and_partial_gas_layer_fits(points):
    one_pass_x_data = []
    run_name_list = []
    points_by_run_name = []
    for point in points:
        if point.theta_r_in_degrees not in one_pass_x_data:
            one_pass_x_data.append(point.theta_r_in_degrees)
        if point.run_name not in run_name_list:
            run_name_list.append(point.run_name)
            points_by_run_name.append([point])
        else:
            index = 0.5  # to error if not overridden
            for i in range(len(run_name_list)):
                run_name = run_name_list[i]
                if point.run_name == run_name:
                    index = i
            points_by_run_name[index].append(point)

    phi_r_in_degrees = points[0].phi_r_in_degrees
    n_0 = points[0].n_0
    polarization = points[0].polarization
    photodiode_solid_angle = points[0].photodiode_solid_angle

    n_gas = 1
    n_liquid = 1.69

    semi_empirical_parameters = semi_empirical_fit.fit_parameters(points)
    TSTR_parameters = TSTR_fit.fit_parameters(points)
    gaussian_parameters = gaussian_fit.fit_parameters(points)
    partial_gas_layer_parameters = partial_gas_layer_fit.fit_parameters_gaussian(points)

    for i in range(len(run_name_list)):

        run_name = run_name_list[i]
        points = points_by_run_name[i]
        theta_i_in_degrees = points[0].theta_i_in_degrees

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]
        max_y = max(y_data)

        semi_empirical_y = semi_empirical_fit.BRIDF_plotter(one_pass_x_data,
                                                            phi_r_in_degrees, theta_i_in_degrees, n_0, polarization,
                                                            photodiode_solid_angle,
                                                            semi_empirical_parameters)
        TSTR_y = TSTR_fit.BRIDF_plotter(one_pass_x_data,
                                        phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, TSTR_parameters)
        gaussian_y = gaussian_fit.BRIDF_plotter(one_pass_x_data, theta_i_in_degrees, gaussian_parameters)
        partial_gas_layer_y = partial_gas_layer_fit.BRIDF_plotter(gaussian_fit.BRIDF_all_parameters, one_pass_x_data,
                                                                  phi_r_in_degrees,
                                                  theta_i_in_degrees, n_gas, n_liquid, polarization,
                                                  photodiode_solid_angle, [partial_gas_layer_parameters[0],
                                                                           partial_gas_layer_parameters[1],
                                                                           partial_gas_layer_parameters[2]],
                                                                  partial_gas_layer_parameters[3])

        plt.figure()
        plt.title(run_name)
        plt.scatter(x_data, y_data, marker="x", color="r", label="experimental")
        plt.plot(one_pass_x_data, semi_empirical_y, label="Semi-Empirical Fit")
        plt.plot(one_pass_x_data, TSTR_y, label="TSTR Fit")
        plt.plot(one_pass_x_data, partial_gas_layer_y, label="Partial Gas Layer Fit")

        plt.plot(one_pass_x_data, gaussian_y, label="Gaussian Fit")

        rho_L_ = semi_empirical_parameters[0]
        n_ = semi_empirical_parameters[1]
        K_ = semi_empirical_parameters[2]
        gamma_ = semi_empirical_parameters[3]
        rho_L = TSTR_parameters[0]
        n = TSTR_parameters[1]
        gamma = TSTR_parameters[2]

        sigma = gaussian_parameters[0]
        R_1 = gaussian_parameters[1]
        R_2 = gaussian_parameters[2]

        string = "theta_i: " + str(theta_i_in_degrees) + "\n\nSemi-Empirical Parameters:\nrho_L: " + str(rho_L_) + \
                 "\nn: " + str(n_) + "\nK: " + str(K_) + "\ngamma: " + str(gamma_) + \
                 "\n\nTSTR Parameters:\nrho_L: " + str(rho_L) + "\nn: " + str(n) + "\ngamma: " + str(gamma) + \
                 "\n\nGaussian Parameters:\nsigma: " + str(sigma) + "\nR_1: " + str(R_1) + "\nR_2: " + str(R_2) + \
                 "\n\nPartial Gas Layer Parameters:\nsigma: " + str(partial_gas_layer_parameters[0]) + "\nR_1: " + \
                 str(partial_gas_layer_parameters[1]) + "\nR_2: " + str(partial_gas_layer_parameters[2]) + "\nx: " + \
                 str(partial_gas_layer_parameters[3])
        axes = plt.gca()
        axes.set_ylim([-0.4 * max_y, 1.2 * max_y])
        plt.legend()
        plt.xlabel("viewing angle (degrees)")
        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.annotate(string, xy=(0.05, 0.55), xycoords='axes fraction', size=6)

    plt.show()


def plot_with_reff_polynomial_fit(points):
    one_pass_x_data = []
    run_name_list = []
    points_by_run_name = []
    for point in points:
        if point.theta_r_in_degrees not in one_pass_x_data:
            one_pass_x_data.append(point.theta_r_in_degrees)
        if point.run_name not in run_name_list:
            run_name_list.append(point.run_name)
            points_by_run_name.append([point])
        else:
            index = 0.5  # to error if not overridden
            for i in range(len(run_name_list)):
                run_name = run_name_list[i]
                if point.run_name == run_name:
                    index = i
            points_by_run_name[index].append(point)

    parameters = reff_polynomial_fit.fit_parameters(points)

    for i in range(len(run_name_list)):

        run_name = run_name_list[i]
        points = points_by_run_name[i]

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]
        max_y = max(y_data)

        y = reff_polynomial_fit.BRIDF_plotter(one_pass_x_data, parameters)

        plt.figure()
        plt.title(run_name)
        plt.scatter(x_data, y_data, marker="x", color="r", label="experimental")
        plt.plot(one_pass_x_data, y, label="REFF polynomial fit")

        c_1 = parameters[0]
        c_2 = parameters[1]
        c_3 = parameters[2]

        string = "theta_i: 0" + "\n\nc_1: " + str(c_1) + "\n\nc_2: " + str(c_2) + "\n\nc_3: " + str(c_3)

        axes = plt.gca()
        axes.set_ylim([-0.4 * max_y, 1.2 * max_y])
        plt.legend()
        plt.xlabel("viewing angle (degrees)")
        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.annotate(string, xy=(0.05, 0.6), xycoords='axes fraction', size=6)

    plt.show()


def plot_with_TSTR_fit(points, title, make_individual_plots=False, log=True):
    one_pass_x_data = []
    run_name_list = []
    points_by_run_name = []
    for point in points:
        if point.theta_r_in_degrees not in one_pass_x_data:
            one_pass_x_data.append(point.theta_r_in_degrees)
        if point.run_name not in run_name_list:
            run_name_list.append(point.run_name)
            points_by_run_name.append([point])
        else:
            index = 0.5  # to error if not overridden
            for i in range(len(run_name_list)):
                run_name = run_name_list[i]
                if point.run_name == run_name:
                    index = i
            points_by_run_name[index].append(point)
    one_pass_x_data = sorted(one_pass_x_data)

    phi_r_in_degrees = points[0].phi_r_in_degrees

    TSTR_parameters = TSTR_fit.fit_parameters(points)

    if make_individual_plots:

        # make a plot for each angle
        for i in range(len(run_name_list)):

            run_name = run_name_list[i]
            points = points_by_run_name[i]
            n_0 = points[0].n_0
            polarization = points[0].polarization
            theta_i_in_degrees = points[0].theta_i_in_degrees

            x_data = [point.theta_r_in_degrees for point in points]
            y_data = [point.intensity for point in points]

            max_y = max(y_data)

            TSTR_y = TSTR_fit.BRIDF_plotter(one_pass_x_data,
                                            phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, TSTR_parameters)

            plt.figure()
            plt.title(run_name)
            plt.scatter(x_data, y_data, color="r", s=2, label="experimental")
            plt.plot(one_pass_x_data, TSTR_y, label="TSTR Fit")

            rho_L = TSTR_parameters[0]
            n = TSTR_parameters[1]
            gamma = TSTR_parameters[2]

            string = "theta_i: " + str(theta_i_in_degrees) + "\n\nTSTR Parameters:\nrho_L: " + \
                     str(rho_L) + "\nn: " + str(n) + "\ngamma: " + str(gamma)

            axes = plt.gca()
            axes.set_ylim([0, 1.2 * max_y])
            if log:
                axes.set_yscale("log", nonposy='clip')
            plt.legend()
            plt.xlabel("viewing angle (degrees)")
            plt.ylabel("intensity (flux/str)/(input flux)")
            plt.annotate(string, xy=(0.05, 0.6), xycoords='axes fraction', size=6)

        plt.figure()

    current_min_y = 1000
    current_max_y = 0
    for i in range(len(run_name_list)):
        run_name = run_name_list[i]
        points = points_by_run_name[i]
        n_0 = points[0].n_0
        polarization = points[0].polarization
        theta_i_in_degrees = points[0].theta_i_in_degrees

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]

        current_min_y = min(current_min_y, min(y_data))
        current_max_y = max(current_max_y, max(y_data))

        TSTR_y = TSTR_fit.BRIDF_plotter(one_pass_x_data,
                                        phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, TSTR_parameters)

        plt.scatter(x_data, y_data, s=2, label=run_name + " experimental")
        plt.plot(one_pass_x_data, TSTR_y, label=run_name + " TSTR Fit")

    plt.title(title)
    axes = plt.gca()
    axes.set_ylim([current_min_y, 1.2 * current_max_y])
    if log:
        axes.set_yscale("log", nonposy='clip')
    plt.legend()
    plt.xlabel("viewing angle (degrees)")
    plt.ylabel("intensity (flux/str)/(input flux)")

    rho_L = TSTR_parameters[0]
    n = TSTR_parameters[1]
    gamma = TSTR_parameters[2]

    string = "TSTR Parameters:\nrho_L: " + \
             str(rho_L) + "\nn: " + str(n) + "\ngamma: " + str(gamma)
    plt.annotate(string, xy=(0.05, 0.6), xycoords='axes fraction', size=6)

    plt.show()


def plot_with_semi_empirical_fit(points, title):
    one_pass_x_data = []
    run_name_list = []
    points_by_run_name = []
    for point in points:
        if point.theta_r_in_degrees not in one_pass_x_data:
            one_pass_x_data.append(point.theta_r_in_degrees)
        if point.run_name not in run_name_list:
            run_name_list.append(point.run_name)
            points_by_run_name.append([point])
        else:
            index = 0.5  # to error if not overridden
            for i in range(len(run_name_list)):
                run_name = run_name_list[i]
                if point.run_name == run_name:
                    index = i
            points_by_run_name[index].append(point)
    one_pass_x_data = sorted(one_pass_x_data)

    phi_r_in_degrees = points[0].phi_r_in_degrees

    semi_empirical_parameters = semi_empirical_fit.fit_parameters(points)

    # make a plot for each angle
    for i in range(len(run_name_list)):

        run_name = run_name_list[i]
        points = points_by_run_name[i]
        n_0 = points[0].n_0
        polarization = points[0].polarization
        theta_i_in_degrees = points[0].theta_i_in_degrees
        photodiode_angle = points[0].photodiode_angular_width

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]

        max_y = max(y_data)

        semi_empirical_y = semi_empirical_fit.BRIDF_plotter(one_pass_x_data,
                                        phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, photodiode_angle, semi_empirical_parameters)

        plt.figure()
        plt.title(run_name)
        plt.scatter(x_data, y_data, color="r", s=2, label="experimental")
        plt.plot(one_pass_x_data, semi_empirical_y, label="Semi-Empirical Fit")

        rho_L = semi_empirical_parameters[0]
        n = semi_empirical_parameters[1]
        K = semi_empirical_parameters[2]
        gamma = semi_empirical_parameters[3]

        string = "theta_i: " + str(theta_i_in_degrees) + "\n\nSemi-Empirical Parameters:\nrho_L: " + \
                 str(rho_L) + "\nn: " + str(n) + "\nK: " + str(K) + "\ngamma: " + str(gamma)

        axes = plt.gca()
        axes.set_ylim([0, 1.2 * max_y])
        plt.legend()
        plt.xlabel("viewing angle (degrees)")
        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.annotate(string, xy=(0.05, 0.6), xycoords='axes fraction', size=6)

    plt.figure()

    current_min_y = 1000
    current_max_y = 0
    for i in range(len(run_name_list)):
        run_name = run_name_list[i]
        points = points_by_run_name[i]
        n_0 = points[0].n_0
        polarization = points[0].polarization
        theta_i_in_degrees = points[0].theta_i_in_degrees

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]

        current_min_y = min(current_min_y, min(y_data))
        current_max_y = max(current_max_y, max(y_data))

        semi_empirical_y = semi_empirical_fit.BRIDF_plotter(one_pass_x_data, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, photodiode_angle, semi_empirical_parameters)

        plt.scatter(x_data, y_data, s=2, label=run_name + " experimental")
        plt.plot(one_pass_x_data, semi_empirical_y, label=run_name + " TSTR Fit")

    plt.title(title)
    axes = plt.gca()
    axes.set_ylim([current_min_y, 1.2 * current_max_y])
    plt.legend()
    plt.xlabel("viewing angle (degrees)")
    plt.ylabel("intensity (flux/str)/(input flux)")

    rho_L = semi_empirical_parameters[0]
    n = semi_empirical_parameters[1]
    K = semi_empirical_parameters[2]
    gamma = semi_empirical_parameters[3]

    string = "\n\nSemi-Empirical Parameters:\nrho_L: " + str(rho_L) + "\nn: " + str(n) + "\nK: " + str(K) + "\ngamma: " + str(gamma)

    plt.annotate(string, xy=(0.05, 0.6), xycoords='axes fraction', size=6)

    plt.show()


def change_theta_i(points, new_theta_i):
    new_points = points[:]
    for point in new_points:
        point.theta_i_in_degrees = new_theta_i
    return new_points


def plot_with_TSTR_fit_and_fitted_angles(points, title):
    new_points = []

    theta_i_list = []
    points_by_run = []
    for point in points:
        if not point.theta_i_in_degrees in theta_i_list:
            theta_i_list.append(point.theta_i_in_degrees)
            points_by_run.append([])
        points_by_run[theta_i_list.index(point.theta_i_in_degrees)].append(point)
    for points in points_by_run:
        intensities = [point.intensity for point in points]
        # determine if there is a sharp peak
        max_point = points[intensities.index(max(intensities))]
        theta_r_peak = max_point.theta_r_in_degrees

        # peak is within 15 degrees of where it's expected
        if np.abs(max_point.theta_r_in_degrees - max_point.theta_i_in_degrees) < 15:
            # use peak as theta_i
            points_1 = change_theta_i(points, theta_r_peak)
            parameters_with_theta_i_peak = TSTR_fit.fit_parameters(points_1)
            theta_r_in_degrees_array = [point.theta_r_in_degrees for point in points]
            point_0 = points_1[0]
            fit_array = TSTR_fit.BRIDF_plotter(theta_r_in_degrees_array, point_0.phi_r_in_degrees, point_0.theta_i_in_degrees,
                                      point_0.n_0, point_0.polarization, parameters_with_theta_i_peak)
            fit_peak_index = fit_array.index(max(fit_array))
            fit_peak = points[fit_peak_index].theta_r_in_degrees
            peak_offset = fit_peak - theta_r_peak
            points = change_theta_i(points, theta_r_peak - peak_offset)
            new_points += points

    plot_with_TSTR_fit(new_points, title)


# uses points to predict a curve
# outputs two plots, one with the points used for fit and the other with the predicted
def TSTR_predict(points_used_for_fit, predicted_points, title):
    return


# returns a list of points from one run of data
# if there are two data points with the same theta_r, it only uses one of them
def make_points(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity_array,
                 wavelength, photodiode_solid_angle, photodiode_angular_width, run_name):
    points = []
    numpoints = len(theta_r_in_degrees_array)
    theta_r_list = [] # used to only include one of each
    for i in range(numpoints):
        if theta_r_in_degrees_array[i] not in theta_r_list:
            theta_r_in_degrees = theta_r_in_degrees_array[i]
            theta_r_list.append(theta_r_in_degrees)
            intensity = intensity_array[i]
            point = Point(theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity,
                     wavelength, photodiode_solid_angle, photodiode_angular_width, run_name)
            points.append(point)
    return points

# returns a list of points from multiple runs of data, finds average and standard deviation
def make_points_std(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity_array,
                 wavelength, photodiode_solid_angle, photodiode_angular_width, run_name):

    theta_r_list = []
    for theta_r in theta_r_in_degrees_array:
        if theta_r not in theta_r_list:
            theta_r_list.append(theta_r)

    # list where the ith element is a list of intensities for the ith theta_r in theta_r_list
    intensities_by_theta_r = []
    for i in theta_r_list:
        intensities_by_theta_r.append([])
    for i in range(len(theta_r_in_degrees_array)):
        pos = theta_r_in_degrees_array.index(theta_r_in_degrees_array[i])
        intensities_by_theta_r[pos].append(intensity_array[i])

    std_array = []
    intensity_averages = []
    for intensity_list in intensities_by_theta_r:
        std_array.append(np.std(intensity_list))
        intensity_averages.append(np.average(intensity_list))

    points = []
    numpoints = len(intensity_averages)
    for i in range(numpoints):
        theta_r_in_degrees = theta_r_list[i]
        intensity = intensity_averages[i]
        std = std_array[i]
        point = Point(theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity, wavelength, photodiode_solid_angle, photodiode_angular_width, run_name, std)
        points.append(point)
    return points


# returns a list of runs which each consist of a list of x_data and a list of y_data
def make_data_by_run(filename, lower_cutoff, upper_cutoff, intensity_factor=1):

    data = np.loadtxt(filename, skiprows=1)
    data_by_run = []

    for i in range(len(data)):
        # if the x_data decreases (new run)
        if i == 0 or data[i][0] < data[i - 1][0]:
            data_by_run.append([[], []])
        if lower_cutoff <= data[i][0] <= upper_cutoff:
            data_by_run[-1][0].append(np.round(data[i][0], 3))
            data_by_run[-1][1].append(intensity_factor * data[i][1])
    return data_by_run


# returns a list of runs which each consist of a list of x_data and a list of y_data
def make_data_by_run_angles_count_down(filename, lower_cutoff, upper_cutoff, intensity_factor=1):

    data = np.loadtxt(filename, skiprows=1)
    data_by_run = []

    for i in range(len(data)):
        # if the x_data decreases (new run)
        if i == 0 or data[i][0] > data[i - 1][0]:
            data_by_run.append([[], []])
        if lower_cutoff <= data[i][0] <= upper_cutoff:
            data_by_run[-1][0].append(np.round(data[i][0], 3))
            data_by_run[-1][1].append(intensity_factor * data[i][1])
    return data_by_run


# returns a list of xdata and a list of ydata
def make_data_all_runs(filename, lower_cutoff, upper_cutoff, intensity_factor=1):
    data = np.loadtxt(filename, skiprows=1)
    x_data = []
    y_data = []

    # this part assumes that if it's less than a fifth of a degree, its an eight of a degree
    if data[1][0] - data[1][1] < 0.2:
        def my_round(x):
            return np.round(x * 8.) / 8.
    elif data[1][0] - data[1][1] < 0.3:
        def my_round(x):
            return np.round(x * 4.) / 4.
    else:
        def my_round(x):
            return np.round(x, 1)

    for i in range(len(data)):
        if lower_cutoff <= data[i][0] <= upper_cutoff:
            x_data.append(my_round(data[i][0]))
            y_data.append(intensity_factor * data[i][1])
    return [x_data, y_data]


def plot_points(points, title, include_individual_plots=False, log=True, show=True, draw_lines=False):
    one_pass_x_data = []
    run_name_list = []
    points_by_run_name = []
    for point in points:
        if point.theta_r_in_degrees not in one_pass_x_data:
            one_pass_x_data.append(point.theta_r_in_degrees)
        if point.run_name not in run_name_list:
            run_name_list.append(point.run_name)
            points_by_run_name.append([point])
        else:
            index = 0.5  # to error if not overridden
            for i in range(len(run_name_list)):
                run_name = run_name_list[i]
                if point.run_name == run_name:
                    index = i
            points_by_run_name[index].append(point)

    if include_individual_plots:
        # make a plot for each angle
        for i in range(len(run_name_list)):
            run_name = run_name_list[i]
            points = points_by_run_name[i]
            theta_i_in_degrees = points[0].theta_i_in_degrees

            x_data = [point.theta_r_in_degrees for point in points]
            y_data = [point.intensity for point in points]

            max_y = max(y_data)

            plt.figure()
            plt.title(run_name)
            plt.scatter(x_data, y_data, color="r", s=2, label="experimental")

            string = "theta_i: " + str(theta_i_in_degrees)

            axes = plt.gca()
            axes.set_ylim([0, 1.2 * max_y])
            plt.legend()
            plt.xlabel("viewing angle (degrees)")
            plt.ylabel("intensity (flux/str)/(input flux)")
            plt.annotate(string, xy=(0.05, 0.6), xycoords='axes fraction', size=6)

    plt.figure()

    min_y = 10**10
    max_y = 0
    for i in range(len(run_name_list)):
        run_name = run_name_list[i]
        points = points_by_run_name[i]

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]

        current_max_y = max(y_data)
        current_min_y = min(y_data)
        max_y = max([max_y, current_max_y])
        min_y = min([min_y, current_min_y])

        plt.scatter(x_data, y_data, s=2, label=run_name)
        if draw_lines:
            plt.plot(x_data, y_data)

    plt.title(title)
    axes = plt.gca()
    axes.set_ylim([0.8 * min_y, 1.2 * max_y])
    if log:
        axes.set_yscale("log", nonposy='clip')
    plt.legend()
    plt.xlabel("viewing angle (degrees)")
    plt.ylabel("intensity (flux/str)/(input flux)")


    if show:
        plt.show()


def plot_points_error_bars(points, title):
    one_pass_x_data = []
    run_name_list = []
    points_by_run_name = []
    for point in points:
        if point.theta_r_in_degrees not in one_pass_x_data:
            one_pass_x_data.append(point.theta_r_in_degrees)
        if point.run_name not in run_name_list:
            run_name_list.append(point.run_name)
            points_by_run_name.append([point])
        else:
            index = 0.5  # to error if not overridden
            for i in range(len(run_name_list)):
                run_name = run_name_list[i]
                if point.run_name == run_name:
                    index = i
            points_by_run_name[index].append(point)

    # make a plot for each angle
    for i in range(len(run_name_list)):
        run_name = run_name_list[i]
        points = points_by_run_name[i]
        theta_i_in_degrees = points[0].theta_i_in_degrees

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]
        err = [point.std for point in points]

        max_y = max(y_data)

        plt.figure()
        plt.title(run_name)
        plt.errorbar(x_data, y_data, yerr=err, linestyle="None", label=run_name)

        string = "theta_i: " + str(theta_i_in_degrees)

        axes = plt.gca()
        axes.set_ylim([0, 1.2 * max_y])
        plt.legend()
        plt.xlabel("viewing angle (degrees)")
        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.annotate(string, xy=(0.05, 0.6), xycoords='axes fraction', size=6)

    max_y = 0
    plt.figure()
    for i in range(len(run_name_list)):
        run_name = run_name_list[i]
        points = points_by_run_name[i]

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]
        err = [point.std for point in points]

        current_max_y = max(y_data)
        max_y = max([max_y, current_max_y])

        plt.title(title)
        plt.errorbar(x_data, y_data, yerr=err, linestyle="None", label=run_name)

        axes = plt.gca()
        axes.set_ylim([0, 1.2 * max_y])
        plt.legend()
        plt.xlabel("viewing angle (degrees)")
        plt.ylabel("intensity (flux/str)/(input flux)")

    plt.show()


def plot_TSTR_fit_error_bars_no_show(points, title):
    color_list = ["r", "g", "b", "m", "c", "y", "k"]

    one_pass_x_data = []
    run_name_list = []
    points_by_run_name = []
    for point in points:
        if point.theta_r_in_degrees not in one_pass_x_data:
            one_pass_x_data.append(point.theta_r_in_degrees)
        if point.run_name not in run_name_list:
            run_name_list.append(point.run_name)
            points_by_run_name.append([point])
        else:
            index = 0.5  # to error if not overridden
            for i in range(len(run_name_list)):
                run_name = run_name_list[i]
                if point.run_name == run_name:
                    index = i
            points_by_run_name[index].append(point)

    max_y = 0
    plt.figure()
    for i in range(len(run_name_list)):
        run_name = run_name_list[i]
        points = points_by_run_name[i]
        n_0 = points[0].n_0
        polarization = points[0].polarization
        theta_i_in_degrees = points[0].theta_i_in_degrees

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]
        err = [point.std for point in points]

        current_max_y = max(y_data)
        max_y = max([max_y, current_max_y])

        plt.title(title)
        plt.errorbar(x_data, y_data, yerr=err, linestyle="None", c=color_list[i])

        fit = TSTR_fit.fit_std(points)

        TSTR_y = TSTR_fit.BRIDF_plotter(one_pass_x_data, 0, theta_i_in_degrees, n_0, polarization, fit)

        plt.plot(one_pass_x_data, TSTR_y, label=run_name, c=color_list[i])

        rho_L = fit[0]
        n = fit[1]
        gamma = fit[2]

        string = str(run_name) + ":\nrho_L: " + str(rho_L) + "\nn: " + str(n) + "\ngamma: " + str(gamma)

        axes = plt.gca()
        axes.set_ylim([0, 1.2 * max_y])
        plt.legend()
        plt.xlabel("viewing angle (degrees)")
        plt.ylabel("intensity (flux/str)/(input flux)")

        plt.annotate(string, xy=(0.05, 0.8 - i / 10.), xycoords='axes fraction', size=6)


def plot_TSTR_fit_one_set_of_parameters(points, parameters, title):
    color_list = ["r", "g", "b", "m", "c", "y", "k", "violet", "darkviolet", "teal"]

    one_pass_x_data = []
    run_name_list = []
    points_by_run_name = []
    for point in points:
        if point.theta_r_in_degrees not in one_pass_x_data:
            one_pass_x_data.append(point.theta_r_in_degrees)
        if point.run_name not in run_name_list:
            run_name_list.append(point.run_name)
            points_by_run_name.append([point])
        else:
            index = 0.5  # to error if not overridden
            for i in range(len(run_name_list)):
                run_name = run_name_list[i]
                if point.run_name == run_name:
                    index = i
            points_by_run_name[index].append(point)

    max_y = 0
    plt.figure()
    for i in range(len(run_name_list)):
        run_name = run_name_list[i]
        points = points_by_run_name[i]

        n_0 = points[0].n_0
        polarization = points[0].polarization
        theta_i_in_degrees = points[0].theta_i_in_degrees

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]
        if points[0].std != None:
            err = [point.std for point in points]
            plt.errorbar(x_data, y_data, yerr=err, linestyle="None", c=color_list[i])
        else:
            plt.scatter(x_data, y_data, s=2, c=color_list[i])

        current_max_y = max(y_data)
        max_y = max([max_y, current_max_y])

        plt.title(title)

        one_pass_x_data = sorted(one_pass_x_data)

        TSTR_y = TSTR_fit.BRIDF_plotter(one_pass_x_data, 0, theta_i_in_degrees, n_0, polarization, parameters)

        plt.plot(one_pass_x_data, TSTR_y, label=run_name, c=color_list[i])

        rho_L = parameters[0]
        n = parameters[1]
        gamma = parameters[2]

        string = "\nrho_L: " + str(rho_L) + "\nn: " + str(n) + "\ngamma: " + str(gamma)

        axes = plt.gca()
        axes.set_ylim([0, 1.2 * max_y])
        axes.set_yscale("log", nonposy='clip')

        plt.legend()
        plt.xlabel("viewing angle (degrees)")
        plt.ylabel("intensity (flux/str)/(input flux)")

        plt.annotate(string, xy=(0.05, 0.4), xycoords='axes fraction', size=10)
    plt.show()


# returns a copy of the points
def copy_points(points):
    copy = []
    for p in points:
        copy.append(Point(p.theta_r_in_degrees, p.phi_r_in_degrees, p.theta_i_in_degrees, p.n_0, p.polarization, p.intensity,
                 p.wavelength, p.photodiode_solid_angle, p.photodiode_angular_width, p.run_name))
    return copy


# returns points with background subtracted
# only returns points which have theta_r the same (when rounded to nearest 8th) for both input lists
# assumes background points have theta_r = 0 at 90 degrees to beam path, 90 in beam path, -90 when blocking beam path
# also adds "with background subtracted" to the runnames
def subtract_background(points, background_points):
    points_copy = copy_points(points)
    # these are relative to the rotation stage angle
    # these have zero in beam path, 180 blocking beam path
    theta_rs_points_rot = [np.round(8 * (180 - p.theta_i_in_degrees - p.theta_r_in_degrees)) / 8. for p in points]
    theta_rs_background = [np.round(8 * (180 - 90 - p.theta_r_in_degrees)) / 8. for p in background_points]

    has_negative = False

    for i in range(len(points)):
        if not theta_rs_points_rot[i] in theta_rs_background:
            print("No background data for theta_rot="+str(theta_rs_points_rot[i]))
        else:
            j = theta_rs_background.index(theta_rs_points_rot[i])
            points_copy[i].intensity = points[i].intensity - background_points[j].intensity
            if points_copy[i].intensity < 0:
                has_negative = True
    for point in points_copy:
        point.run_name = point.run_name + " with background subtracted"
    if has_negative:
        print("Warning: has negative intensities")

    return points_copy


# plots points when x values are decreasing, starts a new line otherwise (removes the annoying slanted lines between plots
def plot_with_restarts(x_data, y_data, label="", color="b"):
    x_data_by_run = []
    y_data_by_run = []
    for i in range(len(x_data)):
        # if we should start a new run
        if i == 0 or x_data[i] < x_data[i - 1]:
            x_data_by_run.append([])
            y_data_by_run.append([])
        x_data_by_run[-1].append(x_data[i])
        y_data_by_run[-1].append(y_data[i])

    for i in range(len(x_data_by_run)):
        if i == 0:
            plt.plot(x_data_by_run[0], y_data_by_run[0], c=color, label=label)
        else:
            plt.plot(x_data_by_run[i], y_data_by_run[i], c=color)
    return
