import radio_beam
from spectral_cube import SpectralCube
from astropy import units as u
import numpy as np

file_in_11 = '../B5_evla_highres/B5_11_VLA_GBT_model_11_22.fits'
file_out_11 = 'data/B5_NH3_11_8arcsec.fits'

file_in_22 = '../B5_evla_highres/B5_22_VLA_GBT_model_11_22.fits'
file_out_22 = 'data/B5_NH3_22_8arcsec.fits'

do_nh3_11 = False
do_nh3_22 = False

do_nh3_22_maps = True

if do_nh3_11:
    cube = SpectralCube.read(file_in_11)
    cube.allow_huge_operations=True
    kcube = cube.to(u.K)
    beam = radio_beam.Beam(major=8*u.arcsec, minor=8*u.arcsec, pa=0*u.deg)
    new_cube = kcube.convolve_to(beam)
    new_cube.write(file_out_11)

if do_nh3_22:
    cube = SpectralCube.read(file_in_22)
    cube.allow_huge_operations=True
    kcube = cube.to(u.K)
    beam = radio_beam.Beam(major=8*u.arcsec, minor=8*u.arcsec, pa=0*u.deg)
    new_cube = kcube.convolve_to(beam)
    new_cube.write(file_out_22)

if do_nh3_22_maps:
    cube_f = SpectralCube.read(file_out_22)
    cube = cube_f.with_spectral_unit(u.km/u.s, velocity_convention='radio')
    spectral_axis = cube.spectral_axis  
    good_channels = (spectral_axis < 10.71*u.km/u.s) & (spectral_axis > 9.71*u.km/u.s)
    bad_channels = ~good_channels
    masked_cube = cube.with_mask(good_channels[:, np.newaxis, np.newaxis])
    noise_cube = cube.with_mask(bad_channels[:, np.newaxis, np.newaxis])
    mom0 = masked_cube.moment( order=0, axis=0)
    rms = noise_cube.std( axis=0)
