from spectral_cube import SpectralCube as SC
import astropy.units as u
import radio_beam

file_in = 'data/B5_argus_gbtidl_gbtpipe_s01-03_cut_Tmb_smoothed5.fits'
cube = SC.read(file_in)
beam = radio_beam.Beam(major=8*u.arcsec, minor=8*u.arcsec, pa=0*u.deg)
new_cube = cube.convolve_to(beam)

file_out = 'data/B5_argus_gbtidl_gbtpipe_s01-03_cut_Tmb_8arcsec.fits'
file_out_grid = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec.fits'
new_cube.write(file_out, overwrite=True)

new_hd = new_cube.header.copy()
new_hd['CDELT1'] = -2.0 / 3600.
new_hd['CDELT2'] = 2.0 / 3600.

from reproject import reproject_interp

new_image, footprint = reproject_interp(new_cube.hdu, new_hd)

from astropy.io import fits
fits.writeto(file_out_grid, new_image, new_hd)
