import numpy as np
import matplotlib.pyplot as plt
from Point import Point
import TSTR_fit
from plotting import plot_with_TSTR_fit, make_points, make_data_by_run, plot_with_TSTR_fit_and_fitted_angles, plot_points

def get_entry_or_default(array, ind, default=0.):
	if ind >= len(array):
		return default
	else:
		return array[ind]

# solid angle of photodiode:
# 5.435 inches from sample
# 9 mm diameter
# in inches
distance_from_sample_to_photodiode = 5.435

# mm / (mm / inch)
photodiode_radius = (9 / 2.0) / 25.4

photodiode_solid_angle = np.pi * np.power(photodiode_radius, 2) / np.power(distance_from_sample_to_photodiode, 2)

photodiode_angular_width = photodiode_radius / distance_from_sample_to_photodiode

# amps per volt during measurements
sensitivity = 100 * 1e-9

wavelength = 405.
pol_equal = 0.5
n_min = 1.47399 # mineral oil
n_gly = 1.47 # glycerol
n_wat = 1.33 # water
n_gly50 = 1.398 # half glycerol, half water
# rearrange to have data by file, to force same num of entries for each data type?

title = "120 Grit Diffuse Reflector, 10 Degree Avg"
file_data=[] # note that each file may have multiple runs
# file_data.append(["glycerol_tests/glycerol 30,45,60,52.5,75.txt",0.00258 * 100e-6])
file_data.append(["glycerol_tests/50_ 30,45,60,52.5,75.txt",0.00252 * 100e-6])

run_data=[] # note that each file may have multiple runs; format is [theta_i, n, polarization, legend, use_in_plots]
# run_data.append([75.,n_gly,pol_equal,"glycerol, 75 degrees", False])
# run_data.append([52.5,n_gly,pol_equal,"glycerol, 52.5 degrees", False])
# run_data.append([60.,n_gly,pol_equal,"glycerol, 60 degrees", True])
# run_data.append([45.,n_gly,pol_equal,"glycerol, 45 degrees", False])
# run_data.append([30.,n_gly,pol_equal,"glycerol, 30 degrees", False])
run_data.append([75.,n_gly50,pol_equal,"50% glycerol, 75 degrees", False])
run_data.append([52.5,n_gly50,pol_equal,"50% glycerol, 52.5 degrees", False])
run_data.append([60.,n_gly50,pol_equal,"50% glycerol, 60 degrees", True])
run_data.append([45.,n_gly50,pol_equal,"50% glycerol, 45 degrees", False])
run_data.append([30.,n_gly50,pol_equal,"30% glycerol, 30 degrees", False])

run_data_transpose = list(map(list, zip(*run_data)))

theta_i_arr = run_data_transpose[0]
n_arr = run_data_transpose[1]
pol_arr = run_data_transpose[2]
leg_arr = run_data_transpose[3] 
use_arr = run_data_transpose[4] # note that data is read out in reverse order of data taking

sub_avg = False
avg_angle = 10.

points_arr = []
ii=0
for entry in file_data:
	filename = entry[0]
	flux_i = entry[1]
	# product of this and measured voltage is (flux/str)/flux_i, flux in units of amps
	# intensity_factor * V = (V * sensitivity / photodiode_solid_angle) / flux_i
	intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)
	data_arr_i = make_data_by_run(filename, -35, 90, intensity_factor)
	for data in data_arr_i:
		if not get_entry_or_default(use_arr, ii, default=False):
			ii += 1
			continue
		if sub_avg:
			y_avg = TSTR_fit.average_by_angle(data[0], data[1], avg_angle)
			y_sub = [data_i-avg_i for data_i, avg_i in zip(data[1],y_avg)]
			data_add = [data[0],y_avg]
		else:
			data_add = data
		leg_i = get_entry_or_default(leg_arr, ii, default="")
		theta_i = get_entry_or_default(theta_i_arr, ii, default=0.)
		n_i = get_entry_or_default(n_arr, ii, default=1.)
		pol_i = get_entry_or_default(pol_arr, ii, default=0.5)
		points_arr += make_points(data_add[0], 0, theta_i, n_i, pol_i, data_add[1], wavelength,
							 photodiode_solid_angle, photodiode_angular_width, leg_i)
		ii += 1

points_x = [point.theta_r_in_degrees for point in points_arr]
points_y = [point.intensity for point in points_arr]

#print(points_x, points_y)
avg_angle = 3.7
fit_params_angle_float = TSTR_fit.fit_parameters_and_angle(points_arr)
fit_params = TSTR_fit.fit_parameters(points_arr, p0=[0.6, 1.3, 0.5], average_angle=avg_angle)
print("Fit parameters: ", fit_params)
#print(fit_params_angle_float)
fit_model = TSTR_fit.BRIDF_plotter(points_x, 0, 60., n_gly50, 0.5, fit_params, average_angle=0.)
fit_model_avg = TSTR_fit.BRIDF_plotter(points_x, 0, 60., n_gly50, 0.5, fit_params, average_angle=avg_angle)
#fit_angle_float = TSTR_fit.BRIDF_plotter(points_x, 0, fit_params_angle_float[0], n_gly50, 0.5, fit_params_angle_float[1:], average_angle=3.7)
model_by_eye = TSTR_fit.BRIDF_plotter(points_x, 0, 60., n_gly50, 0.5, [0.6, 1.27, 0.5], average_angle=avg_angle)
model_by_eye_no_avg = TSTR_fit.BRIDF_plotter(points_x, 0, 60., n_gly50, 0.5, [0.6, 1.27, 0.5], average_angle=0.)
#print(fit)	

res_fit_model = np.sum((np.array(points_y) - np.array(fit_model))**2)
res_fit_model_avg = np.sum((np.array(points_y) - np.array(fit_model_avg))**2)
res_model_by_eye = np.sum((np.array(points_y) - np.array(model_by_eye))**2)
res_model_by_eye_no_avg = np.sum((np.array(points_y) - np.array(model_by_eye_no_avg))**2)
print("residual from fit: ", res_fit_model, ", residual from fit after avg: ", res_fit_model_avg)
print("residual from est by eye: ", res_model_by_eye_no_avg, ", residual by eye after avg: ", res_model_by_eye)

#TSTR_fit.fitter()

plt.scatter(points_x, points_y, s=5, label="data")
plt.semilogy(points_x, fit_model, label="model, fit to data")
#plt.semilogy(points_x, fit_angle_float, label="model, incident angle adj")
plt.semilogy(points_x, model_by_eye, label="model, by eye")
plt.xlabel("viewing angle (degrees)")
plt.ylabel("intensity (flux/str)/(input flux)")
plt.legend()
plt.show()
#plot_points(points_arr, title, log=False, show=False, draw_lines=True)
#plot_points(points_arr, "120 Grit Diffuse Reflector", log=False, show=False)
#plt.show()
