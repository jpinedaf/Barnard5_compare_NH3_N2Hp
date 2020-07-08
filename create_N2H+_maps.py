from spectral_cube import SpectralCube as SC
import astropy.units as u
from astropy.io import fits
import numpy as np
import GAS

from config import file_N2Hp_base_erode, file_N2Hp_base_erode_rms, \
            file_N2Hp_base_erode_TdV, file_N2Hp_base, file_in_N2Hp,\
            file_NH3_11, file_NH3_22, file_NH3_11_match, file_NH3_22_match,\
            file_NH3_11_match_TdV, file_NH3_11_match_rms,\
            file_NH3_22_match_TdV, file_NH3_22_match_rms

# 
# From GAS survey
#
do_baseline_N2Hp = False
if do_baseline_N2Hp:
    GAS.baseline.rebaseline(file_in_N2Hp, blorder=1, 
               baselineRegion=[slice(0, 12, 1), slice(25, 85, 1), 
                               slice(119, 145, 1), slice(177, 199, 1)],
               windowFunction=None, blankBaseline=False,
               flagSpike=False, v0=None, trimEdge=True)


do_trim_N2Hp = False
if do_trim_N2Hp:
    from skimage.morphology import disk,erosion
    cube = SC.read(file_N2Hp_base)
    #
    spectral_axis = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio').spectral_axis  
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
    rms_mask = (rms < 0.5)
    rms_mask &= erosion(rms_mask, disk(5))
    data = cube.unmasked_data[:, :, :] * rms_mask
    data[data[:, :, :] == 0.0] = np.nan
    subcube = SC(data=data, wcs=cube.wcs, header=cube.header)[:, 23:226, 11:166] * u.K
    #
    signal_cube = subcube.with_mask(good_channels[:, np.newaxis, np.newaxis])
    rms_cube = subcube.with_mask(bad_channels[:, np.newaxis, np.newaxis])
    # calcualate rms and integrated itensity maps
    rms = rms_cube.std(axis=0)
    TdV = signal_cube.moment(order=0, axis=0).to(u.K * u.km/u.s)
    # write out the files to be used
    TdV.hdu.writeto(file_N2Hp_base_erode_TdV, overwrite=True)
    rms.hdu.writeto(file_N2Hp_base_erode_rms, overwrite=True)
    subcube.write(file_N2Hp_base_erode, overwrite=True)


do_match_11 = False
if do_match_11:
    cube = SC.read(file_NH3_11)
    cube.allow_huge_operations=True
    # maybe use the same in reproject, to avoid the spectral resampling
    hd_N2Hp = fits.getheader(file_N2Hp_base_erode)
    hd_NH3_11 = cube.header
    key_list = ['NAXIS1', 'NAXIS2', 'BMAJ', 'BMIN', 'BPA', 'CRPIX1', 'CRPIX2',
                'CDELT1', 'CDELT2', 'CUNIT1', 'CUNIT2', 'CTYPE1', 'CTYPE2',
                'CRVAL1', 'CRVAL2']
    for key_i in key_list:
        hd_NH3_11[key_i] = hd_N2Hp[key_i]
    hd_NH3_11.remove('PV2_1', ignore_missing=True)
    hd_NH3_11.remove('PV2_2', ignore_missing=True)
    # print(hd_N2Hp)
    bin_cube = cube.reproject(hd_NH3_11)
    bin_cube.write(file_NH3_11_match, overwrite=True)


do_match_22 = False
if do_match_22:
    cube = SC.read(file_NH3_22)
    # maybe use the same in reproject, to avoid the spectral resampling
    hd_N2Hp = fits.getheader(file_N2Hp_base_erode)
    hd_NH3_22 = cube.header
    key_list = ['NAXIS1', 'NAXIS2', 'BMAJ', 'BMIN', 'BPA', 'CRPIX1', 'CRPIX2',
                'CDELT1', 'CDELT2', 'CUNIT1', 'CUNIT2', 'CTYPE1', 'CTYPE2',
                'CRVAL1', 'CRVAL2']
    for key_i in key_list:
        hd_NH3_22[key_i] = hd_N2Hp[key_i]
    hd_NH3_22.remove('PV2_1', ignore_missing=True)
    hd_NH3_22.remove('PV2_2', ignore_missing=True)
    # print(hd_N2Hp)
    bin_cube = cube.reproject(hd_NH3_22)
    bin_cube.write(file_NH3_22_match, overwrite=True)


do_TdV_11 = False
if do_TdV_11:
    cube_f = SC.read(file_NH3_11_match)
    # cube = cube_f.with_spectral_unit(u.km/u.s, velocity_convention='radio')
    if cube_f.shape[0] > 1018:
        cube = (cube_f[12:-1, :, :]).with_spectral_unit(u.km/u.s,
                                                        velocity_convention='radio')
    else:
        cube = cube_f.with_spectral_unit(u.km/u.s, velocity_convention='radio')
    #
    spectral_axis = cube.spectral_axis  
    good_channels = ((spectral_axis > -9.8*u.km/u.s) & (spectral_axis < -8.6*u.km/u.s)) |\
                    ((spectral_axis > 1.7*u.km/u.s) & (spectral_axis < 3.6*u.km/u.s)) |\
                    ((spectral_axis > 9.3*u.km/u.s) & (spectral_axis < 11.3*u.km/u.s)) |\
                    ((spectral_axis > 17.2*u.km/u.s) & (spectral_axis < 18.8*u.km/u.s)) |\
                    ((spectral_axis > 29.0*u.km/u.s) & (spectral_axis < 30.5*u.km/u.s))
    bad_channels = ~good_channels
    # mask original cube
    signal_cube = cube.with_mask(good_channels[:, np.newaxis, np.newaxis])
    rms_cube = cube.with_mask(bad_channels[:, np.newaxis, np.newaxis])
    # calcualate rms and integrated itensity maps
    rms = rms_cube.std(axis=0)
    TdV = signal_cube.moment(order=0, axis=0).to(u.K * u.km/u.s)
    # write out the files to be used
    TdV.hdu.writeto(file_NH3_11_match_TdV, overwrite=True)
    rms.hdu.writeto(file_NH3_11_match_rms, overwrite=True)


do_TdV_22 = False
if do_TdV_22:
    cube_f = SC.read(file_NH3_22_match)
    cube = cube_f.with_spectral_unit(u.km/u.s, velocity_convention='radio')
    #
    spectral_axis = cube.spectral_axis  
    good_channels = (spectral_axis < 10.71*u.km/u.s) & (spectral_axis > 9.71*u.km/u.s)
    bad_channels = ~good_channels
    # mask original cube
    signal_cube = cube.with_mask(good_channels[:, np.newaxis, np.newaxis])
    rms_cube = cube.with_mask(bad_channels[:, np.newaxis, np.newaxis])
    # calcualate rms and integrated itensity maps
    rms = rms_cube.std(axis=0)
    TdV = signal_cube.moment(order=0, axis=0).to(u.K * u.km/u.s)
    # write out the files to be used
    TdV.hdu.writeto(file_NH3_22_match_TdV, overwrite=True)
    rms.hdu.writeto(file_NH3_22_match_rms, overwrite=True)