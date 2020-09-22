import numpy as np
#import astropy.units as u

file_in_N2Hp = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_full.fits'
file_N2Hp_base = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_full_rebase1.fits'
file_N2Hp_base_erode = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_full_rebase1_erode.fits'
file_N2Hp_base_erode_TdV = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_full_rebase1_erode_TdV.fits'
file_N2Hp_base_erode_rms = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_full_rebase1_erode_rms.fits'

thinFile_N2Hp = 'fit_files/B5_N2Hp_thin_fittedParameters_snr3_full.fits'
thickFile_N2Hp = 'fit_files/B5_N2Hp_thick_fittedParameters_snr3_full.fits'
tpeakFile_N2Hp = 'data/B5_N2Hp_Tpeak_full.fits'
maskFile_N2Hp = 'fit_files/B5_N2Hp_mask_full.fits'
mergedFile_N2Hp = 'fit_files/B5_N2Hp_merged_fittedParameters_full.fits'
mergedFile_N2Hp_Vlsr = 'fit_files/B5_N2Hp_merged_Vlsr_full.fits'
mergedFile_N2Hp_eVlsr = 'fit_files/B5_N2Hp_merged_eVlsr_full.fits'
mergedFile_N2Hp_sigma = 'fit_files/B5_N2Hp_merged_sigma_v_full.fits'
mergedFile_N2Hp_esigma = 'fit_files/B5_N2Hp_merged_esigma_v_full.fits'
mergedFile_N2Hp_Tex = 'fit_files/B5_N2Hp_merged_Tex_full.fits'
mergedFile_N2Hp_eTex = 'fit_files/B5_N2Hp_merged_eTex_full.fits'
mergedFile_N2Hp_tau = 'fit_files/B5_N2Hp_merged_tau_full.fits'
mergedFile_N2Hp_etau = 'fit_files/B5_N2Hp_merged_etau_full.fits'

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
mergedFile_NH3 = 'fit_files/B5_NH3_merged_fittedParameters.fits'
maskFile_NH3 = 'fit_files/B5_NH3_mask.fits'
mergedFile_NH3_Vlsr = 'fit_files/B5_NH3_merged_Vlsr.fits'
mergedFile_NH3_eVlsr = 'fit_files/B5_NH3_merged_eVlsr.fits'
mergedFile_NH3_sigma = 'fit_files/B5_NH3_merged_sigma_v.fits'
mergedFile_NH3_esigma = 'fit_files/B5_NH3_merged_esigma_v.fits'
mergedFile_NH3_Tk = 'fit_files/B5_NH3_merged_Tk.fits'
mergedFile_NH3_eTk = 'fit_files/B5_NH3_merged_eTk.fits'
mergedFile_NH3_Tex = 'fit_files/B5_NH3_merged_Tex.fits'
mergedFile_NH3_eTex = 'fit_files/B5_NH3_merged_eTex.fits'
mergedFile_NH3_N_NH3 = 'fit_files/B5_NH3_merged_N_NH3.fits'
mergedFile_NH3_eN_NH3 = 'fit_files/B5_NH3_merged_eN_NH3.fits'

distance = 300.  # pc
#clev_N2Hp = np.arange(1, 5) * 2
clev_N2Hp = np.arange(5, 25, 5) * 0.4
#clev_NH3_11 = np.arange(1, 5) * 3
clev_NH3_11 = np.array([8, 16, 32, 62]) * 0.3
#clev_NH3_22 = np.arange(1, 5) * 0.2
clev_NH3_22 = np.arange(5, 15, 3) * 0.05

ra_B5IRS1 = (3 + (47 + 41.548/60.)/60.) * 15 # * u.deg
dec_B5IRS1 = (32 + (51 + 43.57/60.)/60.) # * u.deg
ra_B5Cond1 = (3 + (47 + 38.928/60.)/60.) * 15 # * u.deg
dec_B5Cond1 = (32 + (52 + 15.31/60.)/60.) # * u.deg
ra_B5Cond2 = (3 + (47 + 41.627/60.)/60.) * 15 # * u.deg
dec_B5Cond2 = (32 + (51 + 56.81/60.)/60.) # * u.deg
ra_B5Cond3 = (3 + (47 + 42.778/60.)/60.) * 15 # * u.deg
dec_B5Cond3 = (32 + (51 + 30.31/60.)/60.) # * u.deg
