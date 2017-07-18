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


def plot_with_TSTR_fit(points, title):
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
    polarization = points[0].polarization

    TSTR_parameters = TSTR_fit.fit_parameters(points)

    # make a plot for each angle
    for i in range(len(run_name_list)):

        run_name = run_name_list[i]
        points = points_by_run_name[i]
        n_0 = points[0].n_0
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
        plt.legend()
        plt.xlabel("viewing angle (degrees)")
        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.annotate(string, xy=(0.05, 0.6), xycoords='axes fraction', size=6)

    plt.figure()
    for i in range(len(run_name_list)):

        run_name = run_name_list[i]
        points = points_by_run_name[i]
        theta_i_in_degrees = points[0].theta_i_in_degrees

        x_data = [point.theta_r_in_degrees for point in points]
        y_data = [point.intensity for point in points]

        max_y = max(y_data)

        TSTR_y = TSTR_fit.BRIDF_plotter(one_pass_x_data,
                                        phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, TSTR_parameters)

        plt.title(title)
        plt.scatter(x_data, y_data, s=2, label=run_name + " experimental")
        plt.plot(one_pass_x_data, TSTR_y, label=run_name + " TSTR Fit")

        axes = plt.gca()
        axes.set_ylim([0, 1.2 * max_y])
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


def make_points(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity_array,
                 wavelength, photodiode_solid_angle, run_name):
    points = []
    numpoints = len(theta_r_in_degrees_array)
    for i in range(numpoints):
        theta_r_in_degrees = theta_r_in_degrees_array[i]
        intensity = intensity_array[i]
        point = Point(theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity,
                 wavelength, photodiode_solid_angle, run_name)
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
            data_by_run[-1][0].append(np.round(data[i][0], 2))
            data_by_run[-1][1].append(intensity_factor * data[i][1])
    return data_by_run
