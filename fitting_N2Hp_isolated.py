

import astropy.units as u
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from config import file_N2Hp_base_erode, file_N2Hp_base_erode_rms, mergedFile_N2Hp_sigma
from spectral_cube import SpectralCube

import pyspeckit
import matplotlib.pyplot as plt
import numpy as np
plt.ion()

#
# N2Hp (1-0)
#
cube = SpectralCube.read(file_N2Hp_base_erode)

freq_iso = 93176.2604*u.MHz
cube_iso = (cube[0:50, :, :]).with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=freq_iso)
spc = pyspeckit.Cube(cube_iso.hdu)

spc.xarr.convert_to_unit('km/s')
spc.xarr.velocity_convention = 'radio'
spc.xarr.xtype = 'velocity'

#
rms = fits.getdata(file_N2Hp_base_erode_rms)

v_mean = 10.0  # km/s
sigma_v = 0.3  # km/s
Peak = 1.1  # K
my_guess = [Peak, v_mean, sigma_v]

# fit each pixel taking its moment as an initial guess
spc.fiteach(fittype='gaussian',
            guesses=my_guess,  #spc.momentcube,
            errmap=rms,
            signal_cut=3,  # ignore pixels with SNR<3
            blank_value=np.nan,
            start_from_point=(48, 68),
            use_neighbor_as_guess=True,
            multicore=2)
spc.write_fit('fit_files/B5_N2Hp_fit_iso.fits', overwrite=True)

spc.mapplot()
# show the fitted amplitude
#spc.show_fit_param(0, cmap='viridis', vmin=0, vmax=2.0)
spc.show_fit_param(2, cmap='viridis', vmin=0, vmax=0.3)
plt.show()

data, hd = fits.getdata('fit_files/B5_N2Hp_fit_iso.fits', header=True)
key_list = ['NAXIS3', 'CRPIX3', 'CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3']
hd_new = hd.copy()
for key_i in key_list:
    hd_new.remove(key_i)
hd_new['NAXIS'] = 2
hd_new['WCSAXES'] = 2

dv_N2Hp = 0.1 / 2.355  # km/s
dv_iso = np.sqrt(data[2, :, :] ** 2 - dv_N2Hp ** 2)
edv_iso = data[5, :, :]
dv_iso[edv_iso > 0.02] = np.nan
fits.writeto('fit_files/B5_N2Hp_fit_iso_sigma_v.fits', dv_iso, hd_new, overwrite=True)

dv_full = fits.getdata(mergedFile_N2Hp_sigma)
Tp_iso = data[0, :, :]

plt.close()
fig = plt.figure(figsize=(4, 4.5))
ax = fig.add_subplot(111)

im = ax.scatter(dv_full, dv_iso, c=Tp_iso, cmap='YlOrBr', alpha=0.4)
ax.plot([0.05, 0.25], [0.05, 0.25], color='k')
ax.axis('equal')
ax.set_xlim(0.05, 0.25)
ax.set_ylim(0.05, 0.25)
ax.set_xlabel(r'$\sigma_v$ from full spectrum (km s$^{-1}$)')
ax.set_ylabel(r'$\sigma_v$ from isolated component (km s$^{-1}$)')
ax.text(0.23, 0.22, '1:1')
ax.xaxis.set_ticks(np.arange(5e-2, 26e-2, 5e-2))
ax.yaxis.set_ticks(np.arange(5e-2, 26e-2, 5e-2))

#cbaxes = fig.add_axes([0.1, 0.8, 0.8, 0.03])
plt.colorbar(im, orientation='horizontal', label='Line brightness (K)', shrink=0.8)#, cax=cbaxes)
fig.tight_layout()
fig.savefig('figures/compare_sigma_v_N2Hp_iso_full.pdf', dpi=120)
