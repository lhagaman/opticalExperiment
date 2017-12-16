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
# rearrange to have data by file, to force same num of entries for each data type?

title = "120 Grit Diffuse Reflector, 10 Degree Avg"
file_data=[] # note that each file may have multiple runs
#file_data.append(["10_18_diffuse_reflectors/diffuse reflector 1500 grit in mineral oil 30, 45, 60, 75.txt", 0.005830 * 100e-6])
#file_data.append(["10_18_diffuse_reflectors/diffuse reflector 120 grit in mineral oil 30 45 60.txt", 0.005830 * 100e-6])
file_data.append(["11_13_diffuse_reflectors/120 grit 30 45 60 75 11_13.txt",0.00594 * 100e-6])
#file_data.append(["11_13_diffuse_reflectors/120 grit rotated 45 deg 30 45 60 75 11_13.txt",0.00594 * 100e-6])

run_data=[] # note that each file may have multiple runs; format is [theta_i, n, polarization, legend, use_in_plots]
# run_data.append([75.,n_min,pol_equal,"75 degrees", True])
# run_data.append([60.,n_min,pol_equal,"60 degrees", True])
# run_data.append([45.,n_min,pol_equal,"45 degrees", True])
# run_data.append([30.,n_min,pol_equal,"30 degrees", True])
# run_data.append([60.,n_min,pol_equal,"60 degrees, older data", True])
# run_data.append([45.,n_min,pol_equal,"45 degrees, older data", False])
# run_data.append([30.,n_min,pol_equal,"30 degrees, older data", False])
run_data.append([75.,n_min,pol_equal,"75 degrees", True])
run_data.append([60.,n_min,pol_equal,"60 degrees", True])
run_data.append([45.,n_min,pol_equal,"45 degrees", True])
run_data.append([30.,n_min,pol_equal,"30 degrees", True])
# run_data.append([75.,n_min,pol_equal,"75 degrees, sample rotated", False])
# run_data.append([60.,n_min,pol_equal,"60 degrees, sample rotated", True])
# run_data.append([45.,n_min,pol_equal,"45 degrees, sample rotated", False])
# run_data.append([30.,n_min,pol_equal,"30 degrees, sample rotated", False])

run_data_transpose = list(map(list, zip(*run_data)))

theta_i_arr = run_data_transpose[0]
n_arr = run_data_transpose[1]
pol_arr = run_data_transpose[2]
leg_arr = run_data_transpose[3] 
use_arr = run_data_transpose[4] # note that data is read out in reverse order of data taking

sub_avg = True
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
		
plot_points(points_arr, title, log=False, show=False, draw_lines=True)
#plot_points(points_arr, "120 Grit Diffuse Reflector", log=False, show=False)
plt.show()
