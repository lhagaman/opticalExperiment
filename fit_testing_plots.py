# these plot using parameters from the literature, using none of our data.
# This is just to conclude that we've implemented the models correctly.

from plotting import plot_with_semi_empirical_and_gaussian_fits, plot_semi_empirical_components, \
    plot_large_gas_layer_gaussian, plot_partial_gas_layer_gaussian, \
    plot_with_semi_empirical_TSTR_gaussian_and_partial_gas_layer_fits
import numpy as np
import matplotlib.pyplot as plt
from Point import Point
import TSTR_fit
import silva_fit


"""
diffuse = []
lobe = []
spike = []
total = []
angles = [i * 10 for i in range(6)] + [i * 5 for i in range(12, 18)] + [86, 87, 88, 89]
print(angles)

# parameters has the form [rho_L, K, n, gamma], from silvia fit paper
parameters = [0.74, 1.45, 1.7, 0.049]

n_0 = 1
polarization = 0.5
solid_angle = None

for theta_i_in_degrees in angles:
    print("theta_i: ", theta_i_in_degrees)
    diffuse.append(silva_fit.reflectance_diffuse(theta_i_in_degrees * np.pi / 180, n_0, polarization,
                                solid_angle, parameters))
    lobe.append(silva_fit.reflectance_specular_lobe(theta_i_in_degrees * np.pi / 180, n_0, polarization,
                                                 solid_angle, parameters))
    spike.append(silva_fit.reflectance_specular_spike(theta_i_in_degrees * np.pi / 180, n_0, polarization,
                                                 solid_angle, parameters))
    total.append(diffuse[-1] + lobe[-1] + spike[-1])

plt.plot(angles, diffuse, label="diffuse")
plt.plot(angles, lobe, label="specular lobe")
plt.plot(angles, spike, label="specular spike")
plt.plot(angles, total, label="total")

plt.xlabel("angle of incidence (degrees)")
plt.ylabel("reflectance")
plt.legend()
plt.show()
"""

# [rho_L, n, gamma]
# expanded PTFE from paper
parameters = [0.14, 1.56, 0.146]
angles = [i for i in range(-30, 90)]
n_0 = 1
polarization = 0.5

fit_20 = TSTR_fit.BRIDF_plotter(angles, 0, 20, n_0, polarization, parameters)
fit_45 = TSTR_fit.BRIDF_plotter(angles, 0, 45, n_0, polarization, parameters)
fit_65 = TSTR_fit.BRIDF_plotter(angles, 0, 65, n_0, polarization, parameters)

plt.semilogy(angles, fit_20)
plt.semilogy(angles, fit_45)
plt.semilogy(angles, fit_65)

plt.show()





