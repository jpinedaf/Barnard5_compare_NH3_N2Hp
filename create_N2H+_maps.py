from spectral_cube import SpectralCube as SC
import astropy.units as u
from astropy.io import fits
import numpy as np
import GAS

file_in = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec.fits'
file_base = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_rebase1.fits'
file_base_erode = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_rebase1_erode.fits'
file_TdV = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_rebase1_erode_TdV.fits'
file_rms = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_rebase1_erode_rms.fits'

# 
# From GAS survey
#
do_baseline = False
if do_baseline:
    GAS.baseline.rebaseline(file_in, blorder=1, 
               baselineRegion=[slice(0, 12, 1), slice(25, 85, 1), 
                               slice(119, 145, 1), slice(177, 199, 1)],
               windowFunction=None, blankBaseline=False,
               flagSpike=False, v0=None, trimEdge=True)


do_trim = False
if do_trim:
    from skimage.morphology import disk,erosion
    cube = SC.read(file_base)
    #
    spectral_axis = cube.with_spectral_unit(u.km/u.s).spectral_axis  
    good_channels = ((spectral_axis > 1.2*u.km/u.s) & (spectral_axis < 2.5*u.km/u.s)) |\
                    ((spectral_axis > 8.5*u.km/u.s) & (spectral_axis < 11.9*u.km/u.s)) |\
                    ((spectral_axis > 14.5*u.km/u.s) & (spectral_axis < 17.7*u.km/u.s))
    bad_channels = ~good_channels
    # mask original cube
    masked_cube = cube.with_mask(bad_channels[:, np.newaxis, np.newaxis])
    # get rms and mask pixels noisier than 0.5 K
    # after an erosion of the mask
    # save the minimal subcube
    rms = masked_cube.std(axis=0)
    rms_mask = ( rms < 0.5 )
    rms_mask &= erosion( rms_mask, disk(5) )
    data = cube.unmasked_data[:,:,:] * rms_mask
    data[ data[:,:,:] == 0.0 ] = np.nan
    subcube = SC( data=data, wcs=cube.wcs, header=cube.header)[:,23:226,11:166] * u.K
    #
    signal_cube = subcube.with_mask(good_channels[:, np.newaxis, np.newaxis])
    rms_cube = subcube.with_mask(bad_channels[:, np.newaxis, np.newaxis])
    # calcualate rms and integrated itensity maps
    rms = rms_cube.std(axis=0)
    TdV = signal_cube.moment(order=0, axis=0).to(u.K * u.km/u.s)
    # write out the files to be used
    TdV.hdu.writeto(file_TdV, overwrite=True)
    rms.hdu.writeto(file_rms, overwrite=True)
    subcube.write(file_base_erode, overwrite=True)


do_match_22 = True
if do_match_22:
    file_NH3_22 = 'data/B5_NH3_22_8arcsec.fits'
    file_NH3_22_out = 'data/B5_NH3_22_8arcsec_match.fits'
    cube = SC.read(file_NH3_22).with_spectral_unit(u.km/u.s, 
        velocity_convention='radio')
    # maybe use the same in reproject, to avoid the spectral resampling
    hd_N2Hp = fits.getheader(file_base_erode)
    bin_cube = cube.reproject(hd_N2Hp)
    bin_cube.write(file_NH3_22_out)