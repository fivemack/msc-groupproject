from math import exp,pi

arcsec_per_radian = (180*3600)/pi

stefan_boltzmann = 5.67e-8
AU = 1.5e11

asteroid_diameter = 140
telescope_diameter = 0.9
asteroid_temperature = 240
asteroid_emissivity = 0.9

asteroid_range = 1.4*AU

total_asteroid_power = stefan_boltzmann * asteroid_emissivity * 4*pi * (asteroid_diameter/2)**2 * asteroid_temperature**4
print "Power emitted by asteroid is %s W" % total_asteroid_power

power_per_area_at_telescope = total_asteroid_power / (4*pi*asteroid_range**2)
print power_per_area_at_telescope

telescope_aperture_area = pi * (telescope_diameter/2)**2

power_per_telescope = power_per_area_at_telescope * telescope_aperture_area
print "Power from asteroid at telescope is %s W" % power_per_telescope

# statistics to do with the detector
H2RG_pixel_size = 1.8e-5
wavelength = 1e-5
pixels_per_dl = 3
dl_radians = 1.22 * wavelength / telescope_diameter
dl_arcsec = dl_radians * arcsec_per_radian
print "Diffraction limit is %s arcseconds" % dl_arcsec

FOV_pixels = 4096

# we want dl*focal_length=H2RG_pixel_size
focal_length = H2RG_pixel_size / dl_radians
print focal_length
