class Point:
    def __init__(self, theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, intensity,
                 wavelength, photodiode_solid_angle, photodiode_angular_width, run_name, std=None):
        # reflected angle projected to plane of incidence
        self.theta_r_in_degrees = theta_r_in_degrees
        # reflected angle from plane of incidence
        self.phi_r_in_degrees = phi_r_in_degrees
        # incoming angle in plane of incidence
        self.theta_i_in_degrees = theta_i_in_degrees
        # index of refraction outside sample
        self.n_0 = n_0
        # 0 <= polarization <= 1, 0 is s, pendendicular to plane of incidence, 1 is p, parallel to plane of incidence
        self.polarization = polarization
        # measured intensity in units of (flux / str) / flux_i
        self.intensity = intensity
        # wavelength in nanometers
        self.wavelength = wavelength
        # photodiode solid angle, used by semi-empirical fit to determine
        # which points should have the specular spike included
        self.photodiode_solid_angle = photodiode_solid_angle
        # run name to distinguish between different runs and, for example, plot them separately
        self.run_name = run_name
        # this is for fits with a specular spike to see what angles would have the spike included
        self.photodiode_angular_width = photodiode_angular_width
        # this is for fitting when we use average and standard deviation
        self.std = std
