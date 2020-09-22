from spectral_cube import SpectralCube as SC
import astropy.units as u
import radio_beam
from reproject import reproject_interp
from astropy.io import fits
import numpy as np
from astropy.convolution import Gaussian1DKernel

# Load N2H+ cube and then smooth to 8arcsec beam
file_in = 'data/B5_argus_gbtidl_gbtpipe_s01-03_cut_Tmb.fits'
cube = SC.read(file_in)
beam = radio_beam.Beam(major=8*u.arcsec, minor=8*u.arcsec, pa=0*u.deg)
new_cube = cube.convolve_to(beam)
# save file
file_out = 'data/B5_argus_gbtidl_gbtpipe_s01-03_cut_Tmb_8arcsec_full.fits'
new_cube.write(file_out, overwrite=True)

# cube.spectral_axis is np.arange(0,10,0.5) for this example
spectral_axis = new_cube.spectral_axis
fwhm_factor = np.sqrt(8*np.log(2))
current_resolution = spectral_axis[3] - spectral_axis[2]
target_resolution = 0.049 * u.km/u.s
pixel_scale = current_resolution
new_axis = np.arange(spectral_axis[0].value, spectral_axis[-1].value,
                     target_resolution.value) * target_resolution.unit
gaussian_width = ((target_resolution**2 - current_resolution**2)**0.5 /
                  pixel_scale / fwhm_factor)
kernel = Gaussian1DKernel(gaussian_width.value)
smcube = new_cube.spectral_smooth(kernel)
interp_cube = smcube.spectral_interpolate(new_axis,
                                          suppress_smooth_warning=True)

new_hd = interp_cube.header.copy()
new_hd['CDELT1'] = -2.0 / 3600.
new_hd['CDELT2'] = 2.0 / 3600.
# Regrid to smoothed data to match NH3 pixels
new_image, footprint = reproject_interp(interp_cube.hdu, new_hd)
file_out_grid = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_full.fits'
fits.writeto(file_out_grid, new_image, new_hd)
