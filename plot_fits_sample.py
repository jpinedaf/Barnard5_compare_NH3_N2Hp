import pyspeckit
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
import numpy as np
from pyspeckit.spectrum.models import ammonia
from astropy.io import fits
from spectral_cube import SpectralCube as SC
import aplpy

from config import *
#
#import mergedFile_N2Hp, file_N2Hp_base_erode,\
 #   file_N2Hp_base_erode_TdV, file_NH3_11_match_TdV,\
 #   mergedFile_N2Hp, \
 #   mergedFile_NH3, file_NH3_11_match, file_NH3_22_match, \
 #   thinFile_NH3, thickFile_NH3, \
 #   clev_N2Hp, clev_NH3_11, distance

do_N2Hp = True
do_NH3 = True

scale_bar = 10e3 / distance * u.arcsec
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
freqLine11 = 23694.4955 * u.MHz

if do_N2Hp:
    cube = pyspeckit.Cube(file_N2Hp_base_erode)
    cube.xarr.refX = freqLine
    cube.xarr.velocity_convention = 'radio'
    cube.xarr.convert_to_unit('km/s')
    cube.Registry.add_fitter('n2hp_vtau', pyspeckit.models.n2hp.n2hp_vtau_fitter, 4)
    xarr_N2Hp = cube.xarr
    data, hd = fits.getdata(mergedFile_N2Hp, header=True)
    #data[np.isnan(data)] = 0.0
    #fits.writeto('test.fits', data, header=hd, overwrite=True)
    #cube.load_model_fit('test.fits', fittype='n2hp_vtau',
    cube.load_model_fit(mergedFile_N2Hp, fittype='n2hp_vtau',
                        npars=4, npeaks=1, _temp_fit_loc=(x_list[0], y_list[0]))
    cube_sc = (SC.read(file_N2Hp_base_erode)).with_spectral_unit(u.km / u.s,
                                                                velocity_convention='radio')
    model_N2Hp = pyspeckit.models.n2hp.n2hp_vtau_fitter.n_modelfunc

if do_NH3:
    model_NH3 = pyspeckit.spectrum.models.ammonia.cold_ammonia_model().n_modelfunc
    cube11 = pyspeckit.Cube(file_NH3_11_match)
    cube11.xarr.refX = freqLine11
    cube11.xarr.velocity_convention = 'radio'
    cube11.xarr.convert_to_unit('km/s')
    cube11.specfit.Registry.add_fitter('cold_ammonia',
           ammonia.cold_ammonia_model(), 6)
    xarr_NH3 = cube11.xarr
    #thickFile_NH3
    #cube11.load_model_fit(mergedFile_NH3, fittype='cold_ammonia', npars=6, npeaks=1, _temp_fit_loc=(x_list[0], y_list[0]))
    #cube11.load_model_fit('test_fit.fits', fittype='cold_ammonia',
    data2, hd2 = fits.getdata(mergedFile_NH3, header=True)
    #data2[np.isnan(data2)] = 0.0
    #fits.writeto('test2.fits', data2, header=hd2, overwrite=True)
    # data2_hdu = fits.PrimaryHDU(data2, hd2)
    #cube11.load_model_fit('test2.fits', fittype='cold_ammonia',
    cube11.load_model_fit(mergedFile_NH3, fittype='cold_ammonia',
                          npars=6, npeaks=1, _temp_fit_loc=(x_list[0], y_list[0]))
    cube11_sc = (SC.read(file_NH3_11_match)).with_spectral_unit(u.km / u.s,
                                                velocity_convention='radio')
    #cube11_sc[0, :, :] = 0.0

for x_pix, y_pix in zip(x_list, y_list):
    if do_N2Hp:
        cube.plot_spectrum(x_pix, y_pix, plot_fit=True,
                           xmin=-1, xmax=18.5, ypeakscale=1.1)
        plt.savefig('figures/B5_N2Hp_fit_x{0}_y{1}.pdf'.format(x_pix, y_pix))

        fig_a, ax_a = plt.subplots(figsize=(7, 3))
        data_i = cube_sc[:, y_pix, x_pix].value
        model_i = model_N2Hp(data[:, y_pix, x_pix])(xarr_N2Hp)
        ax_a.plot(cube_sc.spectral_axis, data_i,
                  color='k', lw=1, drawstyle='steps-mid')
        ax_a.plot(cube_sc.spectral_axis, model_i,
                  color='red', lw=1)
        ax_a.plot(cube_sc.spectral_axis, model_i - data_i - 0.4,
                  color='gray', lw=1, drawstyle='steps-mid')
        ax_a.set_xlim(0, 18.5)
        ax_a.set_xlabel(r'Velocity (km s$^{-1}$)')
        ax_a.set_ylabel(r'T$_{MB}$ (K)')
        fig_a.savefig('figures/B5_N2Hp_fit_x{0}_y{1}_v2.pdf'.format(x_pix, y_pix))
    # mergedFile_NH3
    if do_NH3:
        cube11.plot_spectrum(x_pix, y_pix, plot_fit=True,
                             xmin=-11, xmax=33, ypeakscale=1.1)
        plt.savefig('figures/B5_NH3_11_fit_x{0}_y{1}.pdf'.format(x_pix, y_pix))

        fig_b, ax_b = plt.subplots(figsize=(7, 3))
        data_i = cube11_sc[:, y_pix, x_pix].value
        model_i = model_NH3(data2[0:6, y_pix, x_pix])(xarr_NH3)
        ax_b.plot(cube11_sc.spectral_axis, data_i,
                  color='k', lw=1, drawstyle='steps-mid')
        ax_b.plot(cube11_sc.spectral_axis, model_i,
                  color='red', lw=1)
        ax_b.plot(cube11_sc.spectral_axis, model_i - data_i - 0.2,
                  color='gray', lw=1, drawstyle='steps-mid')
        ax_b.set_xlim(-12, 33)
        ax_b.set_xlabel(r'Velocity (km s$^{-1}$)')
        ax_b.set_ylabel(r'T$_{MB}$ (K)')
        fig_b.savefig('figures/B5_NH3_fit_x{0}_y{1}_v2.pdf'.format(x_pix, y_pix))


subplot0 = [0.1, 0.1, 0.4, 0.85]
subplot1 = [0.5, 0.1, 0.4, 0.85]

bar_pos0 = [0.4, 0.85, 0.08, 0.025]
bar_pos1 = [0.8, 0.85, 0.08, 0.025]

fig = plt.figure(figsize=(8, 5))
fig0 = aplpy.FITSFigure(file_NH3_11_match_TdV, figure=fig)#,
                        # subplot=subplot0)
fig0.show_colorscale(cmap='Blues', vmin=0, vmax=18)
fig0.set_nan_color(NaN_color)
# Add beam and scalebar
fig0.add_beam(color='black')
fig0.add_scalebar(scale_bar)
fig0.scalebar.set_label(scale_text)
#
# fig1 = aplpy.FITSFigure(file_NH3_11_match_TdV, figure=fig,
#                         subplot=subplot1)
# fig1.show_colorscale(cmap='Reds', vmin=0, vmax=18)
# fig1.set_nan_color(NaN_color)
# # Add beam and scalebar
# fig1.add_beam(color='black')
# fig1.add_scalebar(scale_bar)
# fig1.scalebar.set_label(scale_text)

#
x_world, y_world = fig0.pixel2world(x_list, y_list)
fig0.show_markers(x_world, y_world, color='blue')
# fig0.show_markers(x_world, y_world, color='red')

# No axes labels
fig0.axis_labels.set_xtext('Right Ascension (J2000)')
fig0.axis_labels.set_ytext('Declination (J2000)')
# fig1.axis_labels.hide()
# No tickmarks
# fig1.tick_labels.hide()
# ticks colors
fig0.ticks.set_color('black')
# fig1.ticks.set_color('black')
# contours
# fig0.show_contour(file_N2Hp_base_erode_TdV, levels=clev_N2Hp, colors='k',
#                   linewidths=0.5)
fig0.show_contour(file_NH3_11_match_TdV, levels=clev_NH3_11, colors='k',
                  linewidths=0.5)
# fig0.add_label(0.95, 0.95, r"N$_2$H$^+$(1-0)", relative=True,
#                horizontalalignment='right')
fig0.add_label(0.95, 0.95, r"NH$_3$(1,1)", relative=True,
               horizontalalignment='right')
# fig0.add_colorbar(box=bar_pos0, box_orientation='horizontal',
#                   ticks=[0, 5, 10], axis_label_text=r"(K km s$^{-1}$)")
fig0.add_colorbar(box=bar_pos0, box_orientation='horizontal',
                  ticks=[0, 9, 18], axis_label_text=r"(K km s$^{-1}$)")
fig0.set_system_latex(True)
add_markers_source(fig0, yso_color='yellow', cond_color='0.7')
# fig1.set_system_latex(True)

# fig.savefig('figures/B5_sample_spec.pdf', bbox_inches='tight')
