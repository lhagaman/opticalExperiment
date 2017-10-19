import numpy as np
import matplotlib.pyplot as plt
from TSTR_fit import BRIDF_plotter, BRIDF, BRIDF_radiance_plotter, average_by_angle
from matplotlib.widgets import Slider, Button

log = False
show_averaged = True
angle = 10.
show_all_angles = False

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.5, top=0.95)

theta_r = np.arange(-45, 90, 0.1)
phi_r_0 = 0
theta_i_0 = 45
n_0_0 = 1
polarization_0 = 0.5
rho_l_0 = 0.5
n__0 = 1.5
gamma_0 = 0.05
avg_angle_0 = 4.
s = BRIDF_plotter(theta_r, phi_r_0, theta_i_0, n_0_0, polarization_0, [rho_l_0, n__0, gamma_0])
s1 = average_by_angle(theta_r, s, avg_angle_0)

l, = plt.plot(theta_r, s, lw=2, color='red')
l1, = plt.plot(theta_r, s1, lw=2, color='blue')

if log:
    plt.yscale("log")
    plt.axis([-45., 90., 10e-3, 10e2])
else:
    plt.axis([-45., 90., 0., 10.])

t1 = plt.text(15, 9, "peak angles...")

axcolor = 'lightgoldenrodyellow'

axphi_r = plt.axes([0.25, 0.35, 0.65, 0.03], facecolor=axcolor)
axtheta_i = plt.axes([0.25, 0.3, 0.65, 0.03], facecolor=axcolor)
axn_0 = plt.axes([0.25, 0.25, 0.65, 0.03], facecolor=axcolor)
axpolarization = plt.axes([0.25, 0.2, 0.65, 0.03], facecolor=axcolor)
axrho_l = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
axn = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axgamma = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor=axcolor)
axavg_angle = plt.axes([0.25, 0.00, 0.65, 0.03], facecolor=axcolor)

sphi_r = Slider(axphi_r, "phi_r", -90., 90., valinit=phi_r_0)
stheta_i = Slider(axtheta_i, "theta_i", 0., 90., valinit=theta_i_0)
sn_0 = Slider(axn_0, "n_0", 1., 2., valinit=n_0_0)
spolarization = Slider(axpolarization, "polarization", 0., 1., valinit=polarization_0)
srho_l = Slider(axrho_l, "rho_l", 0., 10., valinit=rho_l_0)
sn = Slider(axn, "n", 1., 2., valinit=n__0)
sgamma = Slider(axgamma, "gamma", 0., 1., valinit=gamma_0)
savg_angle = Slider(axavg_angle,  "average angle", 0, 40, valinit=avg_angle_0)


def update(val):
    phi_r = sphi_r.val
    theta_i = stheta_i.val
    n_0 = sn_0.val
    polarization = spolarization.val
    rho_l = srho_l.val
    n = sn.val
    gamma = sgamma.val
    avg_angle = savg_angle.val
    y_data = BRIDF_plotter(theta_r, phi_r, theta_i, n_0, polarization, [rho_l, n, gamma])
    y1_data = average_by_angle(theta_r, y_data, avg_angle)
    l.set_ydata(y_data)
    l1.set_ydata(y1_data)
    t1.set_text("peak angle is currently " + str(theta_r[y_data.index(max(y_data))]) + "\npeak averaged angle is currently " + str(theta_r[y1_data.index(max(y1_data))]))
    fig.canvas.draw_idle()

sphi_r.on_changed(update)
stheta_i.on_changed(update)
sn_0.on_changed(update)
spolarization.on_changed(update)
srho_l.on_changed(update)
sn.on_changed(update)
sgamma.on_changed(update)
savg_angle.on_changed(update)

resetax = plt.axes([0.8, 0.4, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sphi_r.reset()
    stheta_i.reset()
    sn_0.reset()
    spolarization.reset()
    srho_l.reset()
    sn.reset()
    sgamma.reset()

button.on_clicked(reset)

plt.show()