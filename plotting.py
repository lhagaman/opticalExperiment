import semi_empirical_fit, TSTR_fit, gaussian_fit, large_gas_layer_fit, partial_gas_layer_fit
import numpy as np
import matplotlib.pyplot as plt


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
    n_gas = 1.3
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


def plot_partial_gas_layer_gaussian(points):
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
    n_gas = 1.3
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
        max_y = max(y_data)

        fit_y = partial_gas_layer_fit.BRIDF_plotter(gaussian_fit.BRIDF_all_parameters, one_pass_x_data, phi_r_in_degrees,
                                                  theta_i_in_degrees, n_gas, n_liquid, polarization,
                                                  photodiode_solid_angle, [sigma, R_1, R_2], x)

        plt.figure()
        plt.title(run_name + "\nPartial Gas Layer / Gaussian Fit")
        plt.scatter(x_data, y_data, marker="x", color="r", label="experimental")
        plt.plot(one_pass_x_data, fit_y, label="fit")

        string = "theta_i: " + str(theta_i_in_degrees) + "\nn_liquid: " + str(n_liquid) + "\nn_gas: " + str(n_gas) + \
                 "\n\nx: " + str(x) + "\nsigma: " + str(sigma) + \
                 "\nR_1: " + str(R_1) + "\nR_2: " + str(R_2)

        axes = plt.gca()
        axes.set_ylim([0, 1.2 * max_y])
        plt.legend()
        plt.xlabel("viewing angle (degrees)")
        plt.ylabel("intensity (flux/str)/(input flux)")
        plt.annotate(string, xy=(0.05, 0.7), xycoords='axes fraction', size=6)


