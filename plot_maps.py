import aplpy
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
import numpy as np

import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

from config import file_N2Hp_base_erode, file_N2Hp_base_erode_rms, file_N2Hp_base_erode_TdV,\
            file_NH3_11_match_TdV, file_NH3_22_match_TdV, distance,\
            mergedFile_N2Hp_Vlsr, mergedFile_N2Hp_sigma,\
            mergedFile_NH3_Vlsr, mergedFile_NH3_sigma

scale_bar = 10e3/300 * u.arcsec
scale_text = '10,000 au'

NaN_color='0.9'

do_Tdv = False
do_Vlsr = True

if do_Tdv:
    subplot0=[0.1, 0.1, 0.3, 0.85]
    subplot1=[0.4, 0.1, 0.3, 0.85]
    subplot2=[0.7, 0.1, 0.3, 0.85]

    bar_pos0=[0.3, 0.85, 0.08, 0.025]
    bar_pos1=[0.6, 0.85, 0.08, 0.025]
    bar_pos2=[0.9, 0.85, 0.08, 0.025]

    plt.ion()
    plt.close()

    fig = plt.figure( figsize=(8,3.7))
    fig0=aplpy.FITSFigure(file_N2Hp_base_erode_TdV, figure=fig, subplot=subplot0)

    fig0.show_colorscale(cmap='Blues', vmin=0, vmax=10)
    fig0.set_nan_color(NaN_color)
    # Add beam and scalebar
    fig0.add_beam(color='black')
    fig0.add_scalebar(scale_bar)
    fig0.scalebar.set_label(scale_text)
    #
    fig1=aplpy.FITSFigure(file_NH3_11_match_TdV, figure=fig, subplot=subplot1)

    fig1.show_colorscale(cmap='Reds', vmin=0, vmax=18)
    fig1.set_nan_color(NaN_color)
    # Add beam and scalebar
    fig1.add_beam(color='black')
    fig1.add_scalebar(scale_bar)
    fig1.scalebar.set_label(scale_text)

    fig2=aplpy.FITSFigure(file_NH3_22_match_TdV, figure=fig, subplot=subplot2)

    fig2.show_colorscale(cmap='Reds', vmin=0, vmax=1)
    fig2.set_nan_color(NaN_color)
    # Add beam and scalebar
    fig2.add_beam(color='black')
    fig2.add_scalebar(scale_bar)
    fig2.scalebar.set_label(scale_text)

    # No axes labels
    fig0.axis_labels.set_xtext('Right Ascension (J2000)')
    fig0.axis_labels.set_ytext('Declination (J2000)')
    fig1.axis_labels.hide()
    fig2.axis_labels.hide()
    # No tickmarks
    fig1.tick_labels.hide()
    fig2.tick_labels.hide()
    # ticks colors
    fig0.ticks.set_color('black') 
    fig1.ticks.set_color('black') 
    fig2.ticks.set_color('black') 


    fig0.add_label(0.95, 0.95, r"N$_2$H$^+$(1-0)", relative=True, horizontalalignment='right')
    fig1.add_label(0.95, 0.95, r"NH$_3$(1,1)", relative=True, horizontalalignment='right')
    fig2.add_label(0.95, 0.95, r"NH$_3$(2,2)", relative=True, horizontalalignment='right')

    fig0.add_colorbar(box=bar_pos0, box_orientation='horizontal', ticks=[0, 5, 10],
        axis_label_text=r"(K km s$^{-1}$)")
    fig1.add_colorbar(box=bar_pos1, box_orientation='horizontal', ticks=[0, 9, 18],
        axis_label_text=r"(K km s$^{-1}$)")
    fig2.add_colorbar(box=bar_pos2, box_orientation='horizontal', ticks=[0, 0.5, 1],
        axis_label_text=r"(K km s$^{-1}$)")

    # fig0.colorbar()

    # fig0.colorbar.set_box=bar_pos0, box_orientation='horizontal')
    fig0.set_system_latex(True)
    fig1.set_system_latex(True)
    fig2.set_system_latex(True)

    fig.savefig('figures/B5_TdV_maps.pdf', bbox_inches='tight')

if do_Vlsr:
    subplot0=[0.1, 0.1, 0.3, 0.85]
    subplot1=[0.4, 0.1, 0.3, 0.85]
    subplot2=[0.7, 0.1, 0.3, 0.85]

    bar_pos0=[0.3, 0.85, 0.08, 0.025]
    bar_pos1=[0.6, 0.85, 0.08, 0.025]
    bar_pos2=[0.9, 0.85, 0.08, 0.025]

    plt.ion()
    plt.close()

    fig = plt.figure( figsize=(8,3.7))
    fig0=aplpy.FITSFigure(mergedFile_N2Hp_Vlsr, figure=fig, subplot=subplot0)

    fig0.show_colorscale(cmap='RdYlBu_r', vmin=10, vmax=10.5)
    fig0.set_nan_color(NaN_color)
    # Add beam and scalebar
    fig0.add_beam(color='black')
    fig0.add_scalebar(scale_bar)
    fig0.scalebar.set_label(scale_text)
    #
    fig1=aplpy.FITSFigure(mergedFile_NH3_Vlsr, figure=fig, subplot=subplot1)

    fig1.show_colorscale(cmap='RdYlBu_r', vmin=10, vmax=10.5)
    fig1.set_nan_color(NaN_color)
    # Add beam and scalebar
    fig1.add_beam(color='black')
    fig1.add_scalebar(scale_bar)
    fig1.scalebar.set_label(scale_text)

    # fig2=aplpy.FITSFigure(file_NH3_22_match_TdV, figure=fig, subplot=subplot2)

    # fig2.show_colorscale(cmap='RdYlBu_r', vmin=0, vmax=1)
    # fig2.set_nan_color(NaN_color)
    # # Add beam and scalebar
    # fig2.add_beam(color='black')
    # fig2.add_scalebar(scale_bar)
    # fig2.scalebar.set_label(scale_text)

    # No axes labels
    fig0.axis_labels.set_xtext('Right Ascension (J2000)')
    fig0.axis_labels.set_ytext('Declination (J2000)')
    fig1.axis_labels.hide()
    # fig2.axis_labels.hide()
    # No tickmarks
    fig1.tick_labels.hide()
    # fig2.tick_labels.hide()
    # ticks colors
    fig0.ticks.set_color('black') 
    fig1.ticks.set_color('black') 
    # fig2.ticks.set_color('black') 


    fig0.add_label(0.95, 0.95, r"N$_2$H$^+$(1-0)", relative=True, horizontalalignment='right')
    fig1.add_label(0.95, 0.95, r"NH$_3$(1,1)", relative=True, horizontalalignment='right')
    # fig2.add_label(0.95, 0.95, r"NH$_3$(2,2)", relative=True, horizontalalignment='right')

    fig0.add_colorbar(box=bar_pos0, box_orientation='horizontal', ticks=[0, 5, 10],
        axis_label_text=r"(km s$^{-1}$)")
    fig1.add_colorbar(box=bar_pos1, box_orientation='horizontal', ticks=[0, 9, 18],
        axis_label_text=r"(km s$^{-1}$)")
    # fig2.add_colorbar(box=bar_pos2, box_orientation='horizontal', ticks=[0, 0.5, 1],
    #     axis_label_text=r"(km s$^{-1}$)")

    # fig0.colorbar()

    # fig0.colorbar.set_box=bar_pos0, box_orientation='horizontal')
    fig0.set_system_latex(True)
    fig1.set_system_latex(True)
    # fig2.set_system_latex(True)

    fig.savefig('figures/B5_Vlsr_maps.pdf', bbox_inches='tight')