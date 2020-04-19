from math import exp,pi,sqrt,sin,cos,tan

arcsec_per_radian = (180*3600)/pi
degrees = pi/180 # convert degrees to radians

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
print "Diffraction limit is %.2f arcseconds" % dl_arcsec

FOV_pixels = 4096

# we want dl*focal_length=H2RG_pixel_size
focal_length = pixels_per_dl * H2RG_pixel_size / dl_radians
focal_ratio = focal_length / telescope_diameter
print "Focal length for %s pixels per diffraction-limit is %.2f metres, ratio %.2f" % (pixels_per_dl, focal_length, focal_ratio)

FOV_width_metres = FOV_pixels * H2RG_pixel_size
FOV_width_radians = FOV_pixels * H2RG_pixel_size / focal_length
FOV_width_degrees = FOV_width_radians*arcsec_per_radian/3600

print "Side length of field of view for %s-pixel detector array is %.3f degrees" % (FOV_pixels, FOV_width_degrees)
print "Field of view has %.3g-metre diagonal" % (sqrt(2)*FOV_width_metres)

FOV_area_squaredeg = FOV_width_degrees**2 
sky_area_squaredeg = 360*360/pi

# but we have a 45-degree exclusion window around the Sun
sun_exclusion_degrees = 45

sky_area_squaredeg = sky_area_squaredeg * (1+cos(sun_exclusion_degrees*degrees))/2

fields_per_sky = sky_area_squaredeg / FOV_area_squaredeg

print "%.4g fields-of-view tile the usable sky" % fields_per_sky

