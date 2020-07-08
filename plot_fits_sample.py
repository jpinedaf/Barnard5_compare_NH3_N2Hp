import pyspeckit
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
import numpy as np
from pyspeckit.spectrum.models import ammonia
from astropy.io import fits

from config import mergedFile_N2Hp, file_N2Hp_base_erode,\
    file_N2Hp_base_erode_TdV, file_NH3_11_match_TdV,\
    mergedFile_NH3, file_NH3_11_match, file_NH3_22_match, thickFile_NH3, \
    clev_N2Hp, clev_NH3_11, distance


scale_bar = 10e3/distance * u.arcsec
scale_text = '10,000 au'

NaN_color = '0.9'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
rc('font', **{'family': 'serif'})
rc('text', usetex=True)

plt.ion()

x_list = np.array([63, 62, 45, 49, 48])
y_list = np.array([82, 110, 160, 39, 100])
freqLine = 93173704000.0 * u.Hz

do_N2Hp = True
if do_N2Hp:
    cube = pyspeckit.Cube(file_N2Hp_base_erode)
    cube.xarr.refX = freqLine
    cube.xarr.velocity_convention = 'radio'
    cube.xarr.convert_to_unit('km/s')
    cube.Registry.add_fitter('n2hp_vtau', pyspeckit.models.n2hp.n2hp_vtau_fitter, 4)

do_NH3 = True
if do_NH3:
    #cube.load_model_fit(mergedFile_N2Hp, npars=4, npeaks=1, _temp_fit_loc=(x_list[0], y_list[0]))
    #ff, hd_ff = fits.getdata(file_NH3_11_match, header=True)
    #hd_ff.update(OBJECT='Barnard 5')
    #fits.writeto(file_NH3_11_match, ff, hd_ff, overwrite=True)
    cube11 = pyspeckit.Cube(file_NH3_11_match)
    cube11.specfit.Registry.add_fitter('cold_ammonia',
           ammonia.cold_ammonia_model(), 6)
    #thickFile_NH3
    #cube11.load_model_fit(mergedFile_NH3, fittype='cold_ammonia', npars=6, npeaks=1, _temp_fit_loc=(x_list[0], y_list[0]))
    #cube11.load_model_fit('test_fit.fits', fittype='cold_ammonia',
    cube11.load_model_fit(mergedFile_NH3, fittype='cold_ammonia',
                          npars=6, npeaks=1, _temp_fit_loc=(x_list[0], y_list[0]))

for x_pix, y_pix in zip(x_list, y_list):
    if do_N2Hp:
        cube.plot_spectrum(x_pix, y_pix, plot_fit=True)
        plt.savefig('figures/B5_N2Hp_fit_x{0}_y{1}.pdf'.format(x_pix, y_pix))
    # mergedFile_NH3
    if do_NH3:
        cube11.plot_spectrum(x_pix, y_pix, plot_fit=True)
        plt.savefig('figures/B5_NH3_11_fit_x{0}_y{1}.pdf'.format(x_pix, y_pix))

import aplpy

subplot0 = [0.1, 0.1, 0.4, 0.85]
subplot1 = [0.5, 0.1, 0.4, 0.85]

bar_pos0 = [0.4, 0.85, 0.08, 0.025]
bar_pos1 = [0.8, 0.85, 0.08, 0.025]

fig = plt.figure(figsize=(8, 5))
fig0 = aplpy.FITSFigure(file_N2Hp_base_erode_TdV, figure=fig,
                        subplot=subplot0)
fig0.show_colorscale(cmap='Blues', vmin=0, vmax=10)
fig0.set_nan_color(NaN_color)
# Add beam and scalebar
fig0.add_beam(color='black')
fig0.add_scalebar(scale_bar)
fig0.scalebar.set_label(scale_text)
#
fig1 = aplpy.FITSFigure(file_NH3_11_match_TdV, figure=fig,
                        subplot=subplot1)
fig1.show_colorscale(cmap='Reds', vmin=0, vmax=18)
fig1.set_nan_color(NaN_color)
# Add beam and scalebar
fig1.add_beam(color='black')
fig1.add_scalebar(scale_bar)
fig1.scalebar.set_label(scale_text)

#
x_world, y_world = fig1.pixel2world(x_list, y_list)
fig1.show_markers(x_world, y_world, color='blue')
fig0.show_markers(x_world, y_world, color='red')

# No axes labels
fig0.axis_labels.set_xtext('Right Ascension (J2000)')
fig0.axis_labels.set_ytext('Declination (J2000)')
fig1.axis_labels.hide()
# No tickmarks
fig1.tick_labels.hide()
# ticks colors
fig0.ticks.set_color('black')
fig1.ticks.set_color('black')
# contours
fig0.show_contour(file_N2Hp_base_erode_TdV, levels=clev_N2Hp, colors='k',
                  linewidths=0.5)
fig1.show_contour(file_NH3_11_match_TdV, levels=clev_NH3_11, colors='k',
                  linewidths=0.5)
fig0.add_label(0.95, 0.95, r"N$_2$H$^+$(1-0)", relative=True,
               horizontalalignment='right')
fig1.add_label(0.95, 0.95, r"NH$_3$(1,1)", relative=True,
               horizontalalignment='right')
fig0.add_colorbar(box=bar_pos0, box_orientation='horizontal',
                  ticks=[0, 5, 10], axis_label_text=r"(K km s$^{-1}$)")
fig1.add_colorbar(box=bar_pos1, box_orientation='horizontal',
                  ticks=[0, 9, 18], axis_label_text=r"(K km s$^{-1}$)")
fig0.set_system_latex(True)
fig1.set_system_latex(True)

fig.savefig('figures/B5_sample_spec.pdf', bbox_inches='tight')
