#!/usr/bin/python3
from math import exp,pi,sqrt,sin,cos,tan

arcsec_per_radian = (180*3600)/pi
degrees = pi/180 # convert degrees to radians

stefan_boltzmann = 5.67e-8
planck_h = 6.626e-34
speed_of_light = 2.998e8
AU = 1.5e11

display_intermediate = True

pront = lambda x : x;
if (display_intermediate):
    pront = print

asteroid_diameter = 140
telescope_diameter = 0.9
asteroid_temperature = 240
asteroid_emissivity = 0.9

asteroid_range = 1.4*AU

total_asteroid_power = stefan_boltzmann * asteroid_emissivity * 4*pi * (asteroid_diameter/2)**2 * asteroid_temperature**4
pront("Power emitted by asteroid is %.3g W" % total_asteroid_power )

power_per_area_at_telescope = total_asteroid_power / (4*pi*asteroid_range**2)
telescope_aperture_area = pi * (telescope_diameter/2)**2
bolo_power_per_telescope = power_per_area_at_telescope * telescope_aperture_area
pront( "Bolometric power from asteroid at telescope is %.3g W" % bolo_power_per_telescope )
proportion_in_wavelength_window = 0.107 # this is important and is _not_ the result from planck.py
power_per_telescope = proportion_in_wavelength_window * bolo_power_per_telescope

# statistics to do with the detector
H2RG_pixel_size = 1.8e-5
H2RG_quantum_efficiency = 0.6
H2RG_dark_eps = 200
wavelength = 1e-5
band_short = 8e-6
band_long = 12e-6
photon_energy = planck_h * speed_of_light / wavelength
pixels_per_dl = 3
dl_radians = 1.22 * wavelength / telescope_diameter
dl_arcsec = dl_radians * arcsec_per_radian
pront( "Diffraction limit is %.2f arcseconds" % dl_arcsec )

FOV_pixels = 4096

# we want dl*focal_length=H2RG_pixel_size

focal_length = pixels_per_dl * H2RG_pixel_size / dl_radians
focal_ratio = focal_length / telescope_diameter
pront( "Focal length for %s pixels per diffraction-limit is %.2f metres, ratio %.2f" % (pixels_per_dl, focal_length, focal_ratio) )

FOV_width_metres = FOV_pixels * H2RG_pixel_size
FOV_width_radians = FOV_pixels * H2RG_pixel_size / focal_length
FOV_width_degrees = FOV_width_radians*arcsec_per_radian/3600

pront( "Side length of field of view for %s-pixel detector array is %.3f degrees" % (FOV_pixels, FOV_width_degrees) )
pront( "Field of view has %.3g-metre diagonal" % (sqrt(2)*FOV_width_metres) )

FOV_area_squaredeg = FOV_width_degrees**2 
sky_area_squaredeg = 360*360/pi

# but we have a 45-degree exclusion window around the Sun
sun_exclusion_degrees = 45

sky_area_squaredeg = sky_area_squaredeg * (1+cos(sun_exclusion_degrees*degrees))/2

fields_per_sky = sky_area_squaredeg / FOV_area_squaredeg

pront( "%.4g fields-of-view tile the usable sky" % fields_per_sky )

# the zodiacal light background
zodi_jansky_per_steradian = 3e7 # median figure from [Reach & Morris 2003]
# one Jy is 10^-26 W per square metre per Hz
freq_long = speed_of_light / band_long
freq_short = speed_of_light / band_short
bandwidth_Hz = freq_short - freq_long
zodi_W_m2_sr = zodi_jansky_per_steradian * 1e-26 * bandwidth_Hz
pront( "Zodiacal-light power per square metre aperture per steradian = %.3g W" % zodi_W_m2_sr )
pixel_size_in_sr = (H2RG_pixel_size / focal_length)**2
pront( "Pixel is %.3g sr" % pixel_size_in_sr )
zodi_W_pix = zodi_W_m2_sr * pixel_size_in_sr * telescope_aperture_area
zodi_photons_pix = zodi_W_pix / photon_energy
pront( "Zodiacal light per pixel is %.3g W (%.3g photons/sec)" % (zodi_W_pix, zodi_photons_pix) )

# photons from the source
photons_per_second = power_per_telescope / photon_energy
electrons_per_second = H2RG_quantum_efficiency * photons_per_second

pront( "Faintest source produces %.3g detectable photons (%.3g electrons) per second" % (photons_per_second, electrons_per_second) )

snr = 5
noise = H2RG_dark_eps + zodi_photons_pix*H2RG_quantum_efficiency
pront( "Noise is about %.3g electrons per second" % noise )
exposure = (snr**2*noise) / (electrons_per_second**2)
pront( "Exposure for S/N=5 = %.1f seconds" % exposure )
pront( "Which produces noise=%f  signal=%f  sqrt(noise)=%f" % (noise*exposure, electrons_per_second*exposure, sqrt(noise*exposure)) )

# statistics to do with the spacecraft
pan_time_seconds = 10

sky_scan_time = (exposure + pan_time_seconds) * fields_per_sky

pront( "Scan the sky in %.3g seconds" % sky_scan_time )

