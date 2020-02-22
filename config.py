file_in_N2Hp = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec.fits'
file_N2Hp_base = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_rebase1.fits'
file_N2Hp_base_erode = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_rebase1_erode.fits'
file_N2Hp_base_erode_TdV = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_rebase1_erode_TdV.fits'
file_N2Hp_base_erode_rms = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_rebase1_erode_rms.fits'

thinFile_N2Hp = 'fit_files/B5_N2Hp_thin_fittedParameters_snr3.fits'
thickFile_N2Hp = 'fit_files/B5_N2Hp_thick_fittedParameters_snr3.fits'
tpeakFile_N2Hp = 'data/B5_N2Hp_Tpeak.fits'
maskFile_N2Hp = 'fit_files/B5_N2Hp_mask.fits'
mergedFile_N2Hp = 'fit_files/B5_N2Hp_merged_fittedParameters.fits'
mergedFile_N2Hp_Vlsr = 'fit_files/B5_N2Hp_merged_Vlsr.fits'
mergedFile_N2Hp_eVlsr = 'fit_files/B5_N2Hp_merged_eVlsr.fits'
mergedFile_N2Hp_sigma = 'fit_files/B5_N2Hp_merged_sigma_v.fits'
mergedFile_N2Hp_esigma = 'fit_files/B5_N2Hp_merged_esigma_v.fits'

file_NH3_11 = 'data/B5_NH3_11_8arcsec.fits'
file_NH3_22 = 'data/B5_NH3_22_8arcsec.fits'

file_NH3_11_match = 'data/B5_NH3_11_8arcsec_match.fits'
file_NH3_22_match = 'data/B5_NH3_22_8arcsec_match.fits'
file_NH3_11_match_TdV = 'data/B5_NH3_11_8arcsec_match_TdV.fits'
file_NH3_11_match_Tp = 'data/B5_NH3_11_8arcsec_match_Tp.fits'
file_NH3_11_match_mom1 = 'data/B5_NH3_11_8arcsec_match_mom1.fits'
file_NH3_11_match_rms = 'data/B5_NH3_11_8arcsec_match_rms.fits'
file_NH3_22_match_TdV = 'data/B5_NH3_22_8arcsec_match_TdV.fits'
file_NH3_22_match_rms = 'data/B5_NH3_22_8arcsec_match_rms.fits'

thinFile_NH3 = 'fit_files/B5_NH3_thin_fittedParameters_snr8.fits'
thickFile_NH3 = 'fit_files/B5_NH3_thick_fittedParameters_snr8.fits'
maskFile_NH3 = 'fit_files/B5_NH3_mask.fits'
mergedFile_NH3_Vlsr = 'fit_files/B5_NH3_merged_Vlsr.fits'
mergedFile_NH3_eVlsr = 'fit_files/B5_NH3_merged_eVlsr.fits'
mergedFile_NH3_sigma = 'fit_files/B5_NH3_merged_sigma_v.fits'
mergedFile_NH3_esigma = 'fit_files/B5_NH3_merged_esigma_v.fits'
mergedFile_NH3_Tk = 'fit_files/B5_NH3_merged_Tk.fits'
mergedFile_NH3_eTk = 'fit_files/B5_NH3_merged_eTk.fits'

distance = 300. # pc
import numpy as np
clev_N2Hp = np.arange(1, 5) * 2
clev_NH3_11 = np.arange(1, 5) * 3
clev_NH3_22 = np.arange(1, 5) * 0.2