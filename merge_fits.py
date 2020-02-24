from astropy.io import fits
import numpy as np

from config import thinFile_N2Hp, thickFile_N2Hp, mergedFile_N2Hp,\
    mergedFile_N2Hp_Vlsr, mergedFile_N2Hp_sigma,\
    mergedFile_N2Hp_eVlsr, mergedFile_N2Hp_esigma,\
    mergedFile_NH3_Vlsr, mergedFile_NH3_sigma,\
    mergedFile_NH3_eVlsr, mergedFile_NH3_esigma,\
    mergedFile_NH3_Tk, mergedFile_NH3_eTk,\
    thickFile_NH3

do_N2Hp = True
do_NH3 = True

dv_N2Hp = 0.1 / 2.355 # km/s
dv_NH3 = 0.049 / 2.355  # km/s

if do_N2Hp:
    # PLANE0  = 'TEX     '                                                            
    # PLANE1  = 'TAU     '                                                            
    # PLANE2  = 'CENTER  '                                                            
    # PLANE3  = 'WIDTH   '                                                            
    # PLANE4  = 'eTEX    '                                                            
    # PLANE5  = 'eTAU    '                                                            
    # PLANE6  = 'eCENTER '                                                            
    # PLANE7  = 'eWIDTH  ' 
    thick, hd = fits.getdata(thickFile_N2Hp, header=True)
    thin = fits.getdata(thinFile_N2Hp)
    # mask of tau > 3*sigma_tau
    mask = (thick[ 1, :, :] > 3 * thick[ 5, :, :]) &\
           (thick[ 1, :, :] < 30.)
    merged = mask * thick + (1 - mask) * thin
    # QA
    # Poorly constrained velocity dispersion, errors of 0.02 km/s
    bad = np.broadcast_to(merged[7, :, :] > 0.02, merged.shape)
    merged[bad] = np.nan
    # replace 0 with NaNs
    bad = (merged == 0.)
    if np.sum(bad) != 0:
        merged[bad] = np.nan
    key_list = ['NAXIS3', 'CRPIX3', 'CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3']
    hd_new = hd.copy()
    for key_i in key_list:
        hd_new.remove(key_i)
    hd_new['BUNIT'] = 'km s-1'
    hd_new['NAXIS'] = 2
    hd_new['WCSAXES'] = 2
    # velocity dispersion corrected by channel width
    fits.writeto(mergedFile_N2Hp, merged, hd, overwrite=True)
    fits.writeto(mergedFile_N2Hp_Vlsr, merged[2, :, :], hd_new, overwrite=True)
    fits.writeto(mergedFile_N2Hp_sigma, np.sqrt(merged[3, :, :]**2 - dv_N2Hp**2), hd_new, overwrite=True)
    fits.writeto(mergedFile_N2Hp_eVlsr, merged[6, :, :], hd_new, overwrite=True)
    fits.writeto(mergedFile_N2Hp_esigma, merged[7, :, :] * merged[3, :, :] / np.sqrt(merged[3, :, :]**2 - dv_N2Hp**2), hd_new, overwrite=True)


if do_NH3:
    # PLANE1  = 'TKIN    '                                                            
    # PLANE2  = 'TEX     '                                                            
    # PLANE3  = 'COLUMN  '                                                            
    # PLANE4  = 'SIGMA   '                                                            
    # PLANE5  = 'VELOCITY'                                                            
    # PLANE6  = 'FORTHO  '                                                            
    # PLANE7  = 'eTKIN   '                                                            
    # PLANE8  = 'eTEX    '                                                            
    # PLANE9  = 'eCOLUMN '                                                            
    # PLANE10 = 'eSIGMA  '                                                            
    # PLANE11 = 'eVELOCITY'                                                           
    # PLANE12 = 'eFORTHO '  
    merged, hd = fits.getdata(thickFile_NH3, header=True)
    #thin = fits.getdata(thinFile_NH3)
    # mask of tau > 3*sigma_tau
    # mask = (thick[ 1, :, :] > 3 * thick[ 5, :, :]) &\
    #        (thick[ 1, :, :] < 30.)
    # merged = mask * thick + (1 - mask) * thin
    # QA
    # Poorly constrained velocity dispersion, errors of 0.02 km/s
    bad = np.broadcast_to(merged[9, :, :] > 0.02, merged.shape)
    merged[bad] = np.nan
    # replace 0 with NaNs
    bad = (merged == 0.)
    if np.sum(bad) != 0:
        merged[bad] = np.nan
    key_list = ['NAXIS3', 'CRPIX3', 'CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3']
    hd_new = hd.copy()
    for key_i in key_list:
        hd_new.remove(key_i)
    hd_new['BUNIT'] = 'km s-1'
    hd_new['NAXIS'] = 2
    hd_new['WCSAXES'] = 2
    # fits.writeto(mergedFile_NHp, merged, hd, overwrite=True)
    fits.writeto(mergedFile_NH3_Vlsr, merged[4, :, :], hd_new, overwrite=True)
    fits.writeto(mergedFile_NH3_eVlsr, merged[10, :, :], hd_new, overwrite=True)
    fits.writeto(mergedFile_NH3_sigma, np.sqrt(merged[3, :, :]**2 - dv_NH3**2), hd_new, overwrite=True)
    # fits.writeto(mergedFile_NH3_sigma, merged[3, :, :], hd_new, overwrite=True)
    fits.writeto(mergedFile_NH3_esigma, merged[9, :, :] * merged[3, :, :] / np.sqrt(merged[3, :, :]**2 - dv_NH3**2), hd_new, overwrite=True)
    #
    hd_new['BUNIT'] = 'K'
    fits.writeto(mergedFile_NH3_Tk, merged[0, :, :], hd_new, overwrite=True)
    fits.writeto(mergedFile_NH3_eTk, merged[6, :, :], hd_new, overwrite=True)