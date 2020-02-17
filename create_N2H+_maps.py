from spectral_cube import SpectralCube as SC
import astropy.units as u
from astropy.io import fits
import numpy as np
import GAS

file_in = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec.fits'
# 
# From GAS survey
#
GAS.baseline.rebaseline(file_in, blorder=1, 
               baselineRegion=[slice(0, 12, 1), slice(25, 85, 1), 
                               slice(119, 145, 1), slice(177, 199, 1)],
               windowFunction=None, blankBaseline=False,
               flagSpike=False, v0=None, trimEdge=True)


file_base = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_rebase1.fits'
file_TdV = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_rebase1_TdV.fits'
file_rms = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_rebase1_rms.fits'
cube = SC.read(file_base)

spectral_axis = cube.with_spectral_unit(u.km/u.s).spectral_axis  
good_channels = ((spectral_axis > 1.2*u.km/u.s) & (spectral_axis < 2.5*u.km/u.s)) |
                ((spectral_axis > 8.5*u.km/u.s) & (spectral_axis < 11.9*u.km/u.s)) |
                ((spectral_axis > 14.5*u.km/u.s) & (spectral_axis < 17.7*u.km/u.s))
bad_channels = ~good_channels
masked_cube = cube.with_mask(bad_channels[:, np.newaxis, np.newaxis])
signal_cube = cube.with_mask(good_channels[:, np.newaxis, np.newaxis])
rms = masked_cube.std(axis=0)
TdV = signal_cube.moment(order=0, axis=0)
TdV.hdu.writeto(file_TdV, overwrite=True)
rms.hdu.writeto(file_rms, overwrite=True)
