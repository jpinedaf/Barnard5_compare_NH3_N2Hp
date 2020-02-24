import aplpy
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
import numpy as np
from scipy import stats

import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

from config import file_N2Hp_base_erode, file_N2Hp_base_erode_rms, file_N2Hp_base_erode_TdV,\
            file_NH3_11_match_TdV, file_NH3_22_match_TdV, distance,\
            mergedFile_N2Hp_Vlsr, mergedFile_N2Hp_eVlsr,\
            mergedFile_N2Hp_sigma, mergedFile_N2Hp_esigma,\
            mergedFile_NH3_Vlsr, mergedFile_NH3_eVlsr,\
            mergedFile_NH3_sigma, mergedFile_NH3_esigma,\
            clev_N2Hp, clev_NH3_11, clev_NH3_22

scale_bar = 10e3/300 * u.arcsec
scale_text = '10,000 au'

NaN_color='0.9'
sigma_levels = np.array([0.5, 1.0, 1.5, 2.0])
sigma_levels_l = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
filled_levels = np.hstack([0, sigma_levels]) # 9000 sigmas ~ np.inf
filled_levels_l = np.hstack([0, sigma_levels_l]) # 9000 sigmas ~ np.inf
levels_norm = np.exp(-0.5 * sigma_levels ** 2)[::-1]
levels_norm_f = np.exp(-0.5 * filled_levels ** 2)[::-1]
levels_norm_l = np.exp(-0.5 * filled_levels_l ** 2)[::-1]
color_levels=['#ffffcc', '#ffeda0', '#fed976', '#feb24c', '#fd8d3c', '#fc4e2a', '#e31a1c', '#b10026']
# ['#ece2f0', '#a6bddb', '#1c9099', '#f6eff7','#bdc9e1','#67a9cf','#02818a', '#ece2f0']

do_Tdv = False
do_Vlsr = False
do_sigma_v = False
do_compare_V = False
do_compare_dv = False
do_ratio_dv = True

def my_KDE(data1, data2, xmin=0, xmax=10, ymin=0, ymax=10, weights=None):
    """
    """
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([data1, data2])
    kernel = stats.gaussian_kde(values, weights=weights)
    return X, Y, np.reshape(kernel(positions).T, X.shape)


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
    # contours
    fig0.show_contour(file_N2Hp_base_erode_TdV, levels=clev_N2Hp, colors='k', linewidths=0.5)
    fig1.show_contour(file_NH3_11_match_TdV, levels=clev_NH3_11, colors='k', linewidths=0.5)
    fig2.show_contour(file_NH3_22_match_TdV, levels=clev_NH3_22, colors='k', linewidths=0.5)

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

    fig0.show_contour(file_N2Hp_base_erode_TdV, levels=clev_N2Hp, colors='k', linewidths=0.5)
    fig1.show_contour(file_NH3_11_match_TdV, levels=clev_NH3_11, colors='k', linewidths=0.5)
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

    fig0.add_colorbar(box=bar_pos0, box_orientation='horizontal', ticks=[10, 10.5],
        axis_label_text=r"(km s$^{-1}$)")
    fig1.add_colorbar(box=bar_pos1, box_orientation='horizontal', ticks=[10, 10.5],
        axis_label_text=r"(km s$^{-1}$)")

    fig0.set_system_latex(True)
    fig1.set_system_latex(True)

    fig.savefig('figures/B5_Vlsr_maps.pdf', bbox_inches='tight')


if do_sigma_v:
    subplot0=[0.1, 0.1, 0.3, 0.85]
    subplot1=[0.4, 0.1, 0.3, 0.85]
    subplot2=[0.7, 0.1, 0.3, 0.85]

    bar_pos0=[0.3, 0.85, 0.08, 0.025]
    bar_pos1=[0.6, 0.85, 0.08, 0.025]
    bar_pos2=[0.9, 0.85, 0.08, 0.025]

    plt.ion()
    plt.close()

    fig = plt.figure( figsize=(8,3.7))
    fig0=aplpy.FITSFigure(mergedFile_N2Hp_sigma, figure=fig, subplot=subplot0)

    fig0.show_colorscale(cmap='Blues_r', vmin=0.05, vmax=0.25)
    
    fig0.set_nan_color(NaN_color)
    # Add beam and scalebar
    fig0.add_beam(color='black')
    fig0.add_scalebar(scale_bar)
    fig0.scalebar.set_label(scale_text)
    #
    fig1=aplpy.FITSFigure(mergedFile_NH3_sigma, figure=fig, subplot=subplot1)

    fig1.show_colorscale(cmap='Blues_r', vmin=0.05, vmax=0.25)
    fig1.set_nan_color(NaN_color)
    # Add beam and scalebar
    fig1.add_beam(color='black')
    fig1.add_scalebar(scale_bar)
    fig1.scalebar.set_label(scale_text)

    fig0.show_contour(file_N2Hp_base_erode_TdV, levels=clev_N2Hp, colors='k', linewidths=0.5)
    fig1.show_contour(file_NH3_11_match_TdV, levels=clev_NH3_11, colors='k', linewidths=0.5)
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

    fig0.add_colorbar(box=bar_pos0, box_orientation='horizontal', ticks=[0.05, 0.15, 0.25],
        axis_label_text=r"(km s$^{-1}$)")
    fig1.add_colorbar(box=bar_pos1, box_orientation='horizontal', ticks=[0.05, 0.15, 0.25],
        axis_label_text=r"(km s$^{-1}$)")

    fig0.set_system_latex(True)
    fig1.set_system_latex(True)

    fig.savefig('figures/B5_sigma_v_maps.pdf', bbox_inches='tight')

if do_compare_V:
    plt.ion()
    plt.close()
    fig = plt.figure(figsize=(4,4))
    V_NH3 = fits.getdata(mergedFile_NH3_Vlsr)
    eV_NH3 = fits.getdata(mergedFile_NH3_eVlsr)
    V_N2Hp = fits.getdata(mergedFile_N2Hp_Vlsr)
    eV_N2Hp = fits.getdata(mergedFile_N2Hp_eVlsr)
    gd = np.isfinite(V_NH3 * V_N2Hp)
    plt.errorbar(V_NH3[gd], V_N2Hp[gd], yerr=eV_N2Hp[gd], xerr=eV_NH3[gd], 
        fmt='or', color='red', alpha=0.05)
    xrange = np.array([9.95,10.5])
    plt.plot(xrange, xrange, color='k')
    plt.plot(xrange, xrange+0.05, color='k', ls=':')
    plt.plot(xrange, xrange-0.05, color='k', ls=':')
    plt.xlabel(r"V$_{LSR}$(NH$_3$) (km s$^{-1}$)")
    plt.ylabel(r"V$_{LSR}$(N$_2$H$^+$) (km s$^{-1}$)")
    plt.xlim(xrange)
    plt.ylim(xrange)
    fig.savefig('figures/B5_compare_Vlsr.pdf', bbox_inches='tight')
    # try the KDE
    xx, yy, KDE_vlsr = my_KDE(V_NH3[gd], V_N2Hp[gd], 
                    weights=1./(eV_NH3[gd]*eV_N2Hp[gd]),
                    xmin=xrange[0], xmax=xrange[1], 
                    ymin=xrange[0], ymax=xrange[1])
    fig, ax = plt.subplots(figsize=(4,4))
    levels = levels_norm * KDE_vlsr.max()
    # ax.imshow(np.rot90(KDE_vlsr), cmap='Blues', vmin=0.01, vmax=50, 
    #           extent=[xrange[0], xrange[1], xrange[0], xrange[1]])
    # tkin_kde_colors=['#f6eff7','#bdc9e1','#67a9cf','#02818a']
    # tdust_kde_colors=['#ffffb2','#fecc5c','#fd8d3c','#e31a1c']

    cfset = ax.contourf(xx, yy, KDE_vlsr, colors=color_levels, zorder=2,
                        levels=levels_norm_l * KDE_vlsr.max())
    plt.xlabel(r"V$_{LSR}$(NH$_3$) (km s$^{-1}$)")
    plt.ylabel(r"V$_{LSR}$(N$_2$H$^+$) (km s$^{-1}$)")
    ax.set_xlim(xrange)
    ax.set_ylim(xrange)
    ax.plot(xrange, xrange, color='k', zorder=10)
    ax.plot(xrange, xrange+0.05, color='k', ls=':', zorder=11)
    ax.plot(xrange, xrange-0.05, color='k', ls=':', zorder=12)
    fig.savefig('figures/B5_compare_Vlsr_KDE.pdf', bbox_inches='tight')

if do_compare_dv:
    plt.ion()
    plt.close()
    fig = plt.figure(figsize=(4,4))
    dv_NH3 = fits.getdata(mergedFile_NH3_sigma)
    edv_NH3 = fits.getdata(mergedFile_NH3_esigma)
    dv_N2Hp = fits.getdata(mergedFile_N2Hp_sigma)
    edv_N2Hp = fits.getdata(mergedFile_N2Hp_esigma)
    gd = np.isfinite(dv_NH3 * dv_N2Hp * edv_N2Hp * edv_NH3)
    plt.errorbar(dv_NH3[gd], dv_N2Hp[gd], yerr=edv_N2Hp[gd], xerr=edv_NH3[gd], 
        fmt='or', color='red', alpha=0.01)
    xrange = [0.02, 0.27]
    plt.plot(xrange, xrange)
    plt.xlabel(r"$\sigma_{v}$(NH$_3$) (km s$^{-1}$)")
    plt.ylabel(r"$\sigma_{v}$(N$_2$H$^+$) (km s$^{-1}$)")
    plt.xlim(xrange)
    plt.ylim(xrange)
    fig.savefig('figures/B5_compare_sigma_v.pdf', bbox_inches='tight')
    # try the KDE
    xx, yy, KDE_dv = my_KDE(dv_NH3[gd], dv_N2Hp[gd], 
                    weights=1./(edv_NH3[gd]*edv_N2Hp[gd]),
                    xmin=xrange[0], xmax=xrange[1], 
                    ymin=xrange[0], ymax=xrange[1])
    fig, ax = plt.subplots(figsize=(4,4))
    levels = levels_norm * KDE_dv.max()
    cfset = ax.contourf(xx, yy, KDE_dv, colors=color_levels, zorder=2,
                        levels=levels_norm_l * KDE_dv.max())
    plt.xlabel(r"$\sigma_{v}$(NH$_3$) (km s$^{-1}$)")
    plt.ylabel(r"$\sigma_{v}$(N$_2$H$^+$) (km s$^{-1}$)")
    ax.set_xlim(xrange)
    ax.set_ylim(xrange)
    ax.plot(xrange, xrange, color='k', zorder=10)
    fig.savefig('figures/B5_compare_sigma_v_KDE.pdf', bbox_inches='tight')

if do_ratio_dv:
    plt.ion()
    plt.close()
    dv_NH3, hd = fits.getdata(mergedFile_NH3_sigma, header=True)
    edv_NH3 = fits.getdata(mergedFile_NH3_esigma)
    dv_N2Hp = fits.getdata(mergedFile_N2Hp_sigma)
    edv_N2Hp = fits.getdata(mergedFile_N2Hp_esigma)
    mask_dv = np.isfinite(dv_NH3 * dv_N2Hp * edv_N2Hp * edv_NH3)
    ratio_dv = dv_NH3 / dv_N2Hp
    ratio_dv[~mask_dv] = np.nan
    fig = plt.figure(figsize=(4,3.7))
    fig0=aplpy.FITSFigure(fits.PrimaryHDU(ratio_dv, header=hd), figure=fig)#, subplot=subplot1)
    fig0.show_colorscale(cmap='RdYlBu', vmin=0.5, vmax=1.5)
    fig0.set_nan_color(NaN_color)
    fig0.show_contour(file_NH3_11_match_TdV, levels=clev_NH3_11, colors='k', linewidths=0.5)
    fig0.axis_labels.set_xtext('Right Ascension (J2000)')
    fig0.axis_labels.set_ytext('Declination (J2000)')
    fig0.ticks.set_color('black')
    fig0.add_beam(color='black')
    fig0.add_scalebar(scale_bar)
    fig0.scalebar.set_label(scale_text)
    fig0.add_colorbar(box=[0.50,0.8,0.25,0.05], box_orientation='horizontal', ticks=[0.5, 1, 1.5],
        axis_label_text=r"$\sigma_v$(N$_2$H$^+$)/$\sigma_v$(NH$_3$)")
    fig0.set_system_latex(True)
    fig.savefig('figures/B5_ratio_sigma_v_maps.pdf', bbox_inches='tight')