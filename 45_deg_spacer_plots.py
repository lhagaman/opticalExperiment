import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("45 deg debugging tests.txt", skiprows=1)
data_x = [np.round(point[0], 2) for point in data]
data_y = [point[1] for point in data]

cutoff_data_x = []
cutoff_data_y = []
for i in range(len(data_x)):
    if data_x[i] > -21:
        cutoff_data_x.append(data_x[i])
        cutoff_data_y.append(data_y[i])

data_x_by_run = []
data_y_by_run = []
j = -1
for i in range(len(cutoff_data_x)):
    if cutoff_data_x[i] == -20.5:
        j += 1
        data_x_by_run.append([])
        data_y_by_run.append([])
    data_x_by_run[j].append(cutoff_data_x[i])
    data_y_by_run[j].append(cutoff_data_y[i])

plt.plot(data_x_by_run[0], data_y_by_run[0], color="b", label="with spacer")
plt.plot(data_x_by_run[1], data_y_by_run[1], color="darkturquoise", label="with spacer after changing angles")
plt.plot(data_x_by_run[2], data_y_by_run[2], color="teal", label="with spacer after changing angles again")
plt.plot(data_x_by_run[3], data_y_by_run[3], color="paleturquoise", label="with spacer after replacing sample")
plt.plot(data_x_by_run[4], data_y_by_run[4], color="r", label="without spacer")
plt.plot(data_x_by_run[5], data_y_by_run[5], color="darkorange", label="without spacer after replacing sample")
plt.legend()
plt.xlabel("viewing angle (degrees)")
plt.ylabel("intensity (Volts)")
plt.title("45 degrees in water")
plt.show()
