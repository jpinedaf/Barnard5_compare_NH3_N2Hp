import aplpy
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
import numpy as np

file_TdV = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_rebase1_erode_TdV.fits'
file_rms = 'data/B5_N2H+_1-0_ARGUS_Tmb_8arcsec_rebase1_erode_rms.fits'

distance = 300. # pc
scale_bar = 10e3/300 * u.arcsec
scale_text = '10,000 au'

plt.ion()

fig = plt.figure( figsize=(8,6))

fig0=aplpy.FITSFigure(file_TdV, figure=fig)

fig0.show_colorscale(cmap='Blues')
fig0.set_nan_color('0.99')
# Add beam and scalebar
fig0.add_beam()
fig0.add_scalebar(scale_bar)
fig0.scalebar.set_label(scale_text)