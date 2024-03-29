import aplpy
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
import numpy as np
from scipy import stats
import matplotlib as mpl
from matplotlib import rc
from astropy.constants import k_B
from config import file_N2Hp_base_erode, file_N2Hp_base_erode_rms, file_N2Hp_base_erode_TdV,\
            file_NH3_11_match_TdV, file_NH3_22_match_TdV, distance,\
            mergedFile_N2Hp_Vlsr, mergedFile_N2Hp_eVlsr,\
            mergedFile_N2Hp_sigma, mergedFile_N2Hp_esigma, \
            mergedFile_N2Hp_Tex, mergedFile_N2Hp_eTex, \
            mergedFile_N2Hp_tau, mergedFile_N2Hp_etau,\
            mergedFile_NH3_Vlsr, mergedFile_NH3_eVlsr,\
            mergedFile_NH3_sigma, mergedFile_NH3_esigma, \
            mergedFile_NH3_Tk, mergedFile_NH3_eTk, \
            mergedFile_NH3_Tex, mergedFile_NH3_eTex,\
            clev_N2Hp, clev_NH3_11, clev_NH3_22,\
            ra_B5IRS1, dec_B5IRS1, ra_B5Cond1, dec_B5Cond1,\
            ra_B5Cond2, dec_B5Cond2, ra_B5Cond3, dec_B5Cond3
plt.ion()
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
rc('font', **{'family': 'serif'})
rc('text', usetex=True)

scale_bar = 10e3/distance * u.arcsec
scale_text = '10,000 au'

NaN_color = '0.9'
sigma_levels = np.array([0.5, 1.0, 1.5, 2.0])
sigma_levels_l = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
filled_levels = np.hstack([0, sigma_levels])  # 9000 sigmas ~ np.inf
filled_levels_l = np.hstack([0, sigma_levels_l])  # 9000 sigmas ~ np.inf
levels_norm = np.exp(-0.5 * sigma_levels ** 2)[::-1]
levels_norm_f = np.exp(-0.5 * filled_levels ** 2)[::-1]
levels_norm_l = np.exp(-0.5 * filled_levels_l ** 2)[::-1]
color_levels = ['#ffffcc', '#ffeda0', '#fed976', '#feb24c',
                '#fd8d3c', '#fc4e2a', '#e31a1c', '#b10026']
color_levels2 = ['#deebf7', '#c6dbef', '#9ecae1', '#6baed6',
                 '#4292c6', '#2171b5', '#08519c', '#08306b']
color_N2Hp = '#377eb8'
color_NH3 = '#e41a1c'

do_Tdv = False
do_Vlsr = False
do_Vlsr_diff = False
do_sigma_v = False
do_compare_V = False
do_compare_dv = False
do_ratio_dv = False
do_compare_dv_nt = True


def my_KDE(data1, data2, xmin=0, xmax=10, ymin=0, ymax=10, weights=None):
    """
    """
    x, y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions_kde = np.vstack([x.ravel(), y.ravel()])
    values = np.vstack([data1, data2])
    kernel = stats.gaussian_kde(values, weights=weights)
    return x, y, np.reshape(kernel(positions_kde).T, x.shape)


def sigma_thermal(mu_mol, tk=10*u.K):
    """
    Returns the sound speed for temperature Tk and molecular weight mu.
    This is also used to determine the thermal velocity dispersion of
    a molecular transition.

    """
    return np.sqrt(k_B * tk/(mu_mol * u.u)).to(u.km/u.s)


def add_markers_source(figure, yso_color='yellow', cond_color='black'):
    """

    """
    figure.show_markers(ra_B5IRS1, dec_B5IRS1, marker='*', s=40, layer='YSO',
                        facecolor=yso_color, edgecolor='k', zorder=31)
    figure.show_markers(np.array([ra_B5Cond1, ra_B5Cond2, ra_B5Cond3]),
                        np.array([dec_B5Cond1, dec_B5Cond2, dec_B5Cond3]),
                        marker='o', s=10, layer='Condensations', alpha=0.5,
                        facecolor=cond_color, edgecolor='k', zorder=32)
    return

if do_Tdv:
    subplot0 = [0.1, 0.1, 0.3, 0.85]
    subplot1 = [0.4, 0.1, 0.3, 0.85]
    subplot2 = [0.7, 0.1, 0.3, 0.85]

    bar_pos0 = [0.3, 0.85, 0.08, 0.025]
    bar_pos1 = [0.6, 0.85, 0.08, 0.025]
    bar_pos2 = [0.9, 0.85, 0.08, 0.025]

    plt.ion()
    plt.close()

    fig = plt.figure(figsize=(8, 3.7))
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

    fig2 = aplpy.FITSFigure(file_NH3_22_match_TdV, figure=fig,
                            subplot=subplot2)
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
    fig0.show_contour(file_N2Hp_base_erode_TdV, levels=clev_N2Hp, colors='k',
                      linewidths=0.5)
    fig1.show_contour(file_NH3_11_match_TdV, levels=clev_NH3_11, colors='k',
                      linewidths=0.5)
    fig2.show_contour(file_NH3_22_match_TdV, levels=clev_NH3_22, colors='k',
                      linewidths=0.5)
    fig0.add_label(0.95, 0.95, r"N$_2$H$^+$(1-0)", relative=True,
                   horizontalalignment='right')
    fig1.add_label(0.95, 0.95, r"NH$_3$(1,1)", relative=True,
                   horizontalalignment='right')
    fig2.add_label(0.95, 0.95, r"NH$_3$(2,2)", relative=True,
                   horizontalalignment='right')
    fig0.add_colorbar(box=bar_pos0, box_orientation='horizontal',
                      ticks=[0, 5, 10], axis_label_text=r"(K km s$^{-1}$)")
    fig1.add_colorbar(box=bar_pos1, box_orientation='horizontal',
                      ticks=[0, 9, 18], axis_label_text=r"(K km s$^{-1}$)")
    fig2.add_colorbar(box=bar_pos2, box_orientation='horizontal',
                      ticks=[0, 0.5, 1], axis_label_text=r"(K km s$^{-1}$)")
    fig0.set_system_latex(True)
    fig1.set_system_latex(True)
    fig2.set_system_latex(True)

    fig.savefig('figures/B5_TdV_maps.pdf', bbox_inches='tight')
    add_markers_source(fig0, yso_color='yellow', cond_color='0.7')
    add_markers_source(fig1, yso_color='yellow', cond_color='0.7')
    add_markers_source(fig2, yso_color='yellow', cond_color='0.7')
    fig.savefig('figures/B5_TdV_maps_markers.pdf', bbox_inches='tight')

if do_Vlsr:
    subplot0 = [0.1, 0.1, 0.3, 0.85]
    subplot1 = [0.4, 0.1, 0.3, 0.85]
    subplot2 = [0.7, 0.1, 0.3, 0.85]

    bar_pos0 = [0.3, 0.85, 0.08, 0.025]
    bar_pos1 = [0.6, 0.85, 0.08, 0.025]
    bar_pos2 = [0.9, 0.85, 0.08, 0.025]

    plt.ion()
    plt.close()

    fig = plt.figure(figsize=(8, 3.7))
    fig0 = aplpy.FITSFigure(mergedFile_N2Hp_Vlsr, figure=fig, subplot=subplot0)
    fig0.show_colorscale(cmap='RdYlBu_r', vmin=10, vmax=10.5)
    fig0.set_nan_color(NaN_color)
    # Add beam and scalebar
    fig0.add_beam(color='black')
    fig0.add_scalebar(scale_bar)
    fig0.scalebar.set_label(scale_text)
    #
    fig1 = aplpy.FITSFigure(mergedFile_NH3_Vlsr, figure=fig, subplot=subplot1)
    fig1.show_colorscale(cmap='RdYlBu_r', vmin=10, vmax=10.5)
    fig1.set_nan_color(NaN_color)
    # Add beam and scalebar
    fig1.add_beam(color='black')
    fig1.add_scalebar(scale_bar)
    fig1.scalebar.set_label(scale_text)
    fig0.show_contour(file_N2Hp_base_erode_TdV, levels=clev_N2Hp, colors='k',
                      linewidths=0.5)
    fig1.show_contour(file_NH3_11_match_TdV, levels=clev_NH3_11, colors='k',
                      linewidths=0.5)
    # No axes labels
    fig0.axis_labels.set_xtext('Right Ascension (J2000)')
    fig0.axis_labels.set_ytext('Declination (J2000)')
    fig1.axis_labels.hide()
    # No tickmarks
    fig1.tick_labels.hide()
    # ticks colors
    fig0.ticks.set_color('black') 
    fig1.ticks.set_color('black') 
    #
    fig0.add_label(0.95, 0.95, r"N$_2$H$^+$(1-0)", relative=True,
                   horizontalalignment='right')
    fig1.add_label(0.95, 0.95, r"NH$_3$(1,1)", relative=True,
                   horizontalalignment='right')
    fig0.add_colorbar(box=bar_pos0, box_orientation='horizontal',
                      ticks=[10, 10.5], axis_label_text=r"(km s$^{-1}$)")
    fig1.add_colorbar(box=bar_pos1, box_orientation='horizontal',
                      ticks=[10, 10.5], axis_label_text=r"(km s$^{-1}$)")
    fig0.set_system_latex(True)
    fig1.set_system_latex(True)
    add_markers_source(fig0, yso_color='yellow', cond_color='0.7')
    add_markers_source(fig1, yso_color='yellow', cond_color='0.7')
    fig.savefig('figures/B5_Vlsr_maps.pdf', bbox_inches='tight')

if do_Vlsr_diff:
    bar_pos0 = [0.60, 0.5, 0.20, 0.025]

    plt.ion()
    plt.close()
    v_N2Hp, hd = fits.getdata(mergedFile_N2Hp_Vlsr, header=True)
    v_NH3 = fits.getdata(mergedFile_NH3_Vlsr)
    v_diff = v_N2Hp - v_NH3
    hdu_diff = fits.PrimaryHDU(v_diff, header=hd)
    fig = plt.figure(figsize=(4, 5))
    fig0 = aplpy.FITSFigure(hdu_diff, figure=fig)
    fig0.show_colorscale(cmap='RdBu_r', vmin=-0.1, vmax=0.1)
    fig0.set_nan_color(NaN_color)
    # Add beam and scalebar
    fig0.add_beam(color='black')
    fig0.add_scalebar(scale_bar)
    fig0.scalebar.set_label(scale_text)
    #
    fig0.show_contour(file_NH3_11_match_TdV, levels=clev_NH3_11, colors='k',
                      linewidths=0.5)
    # No axes labels
    fig0.axis_labels.set_xtext('Right Ascension (J2000)')
    fig0.axis_labels.set_ytext('Declination (J2000)')
    # ticks colors
    fig0.ticks.set_color('black')
    add_markers_source(fig0, yso_color='yellow', cond_color='0.7')
    #
    fig0.add_label(0.95, 0.81, r"$\delta V_{\rm LSR}$= ion - neutral", relative=True,
                   horizontalalignment='right')
    fig0.add_colorbar(box=bar_pos0, box_orientation='horizontal',
                      ticks=[-0.1, 0, 0.1], axis_label_text=r"(km s$^{-1}$)")
    fig0.set_system_latex(True)
    fig.savefig('figures/B5_Vlsr_diff_maps.pdf', bbox_inches='tight', dpi=120)

if do_sigma_v:
    subplot0 = [0.1, 0.1, 0.3, 0.85]
    subplot1 = [0.4, 0.1, 0.3, 0.85]
    subplot2 = [0.7, 0.1, 0.3, 0.85]

    bar_pos0 = [0.3, 0.85, 0.08, 0.025]
    bar_pos1 = [0.6, 0.85, 0.08, 0.025]
    bar_pos2 = [0.9, 0.85, 0.08, 0.025]

    plt.ion()
    plt.close()

    fig = plt.figure(figsize=(8, 3.7))
    fig0 = aplpy.FITSFigure(mergedFile_N2Hp_sigma, figure=fig,
                            subplot=subplot0)
    fig0.show_colorscale(cmap='Blues_r', vmin=0.05, vmax=0.25)
    fig0.set_nan_color(NaN_color)
    # Add beam and scalebar
    fig0.add_beam(color='black')
    fig0.add_scalebar(scale_bar)
    fig0.scalebar.set_label(scale_text)
    #
    fig1 = aplpy.FITSFigure(mergedFile_NH3_sigma, figure=fig,
                            subplot=subplot1)
    fig1.show_colorscale(cmap='Blues_r', vmin=0.05, vmax=0.25)
    fig1.set_nan_color(NaN_color)
    # Add beam and scalebar
    fig1.add_beam(color='black')
    fig1.add_scalebar(scale_bar)
    fig1.scalebar.set_label(scale_text)

    fig0.show_contour(file_N2Hp_base_erode_TdV, levels=clev_N2Hp, colors='k',
                      linewidths=0.5)
    fig1.show_contour(file_NH3_11_match_TdV, levels=clev_NH3_11, colors='k',
                      linewidths=0.5)
    # No axes labels
    fig0.axis_labels.set_xtext('Right Ascension (J2000)')
    fig0.axis_labels.set_ytext('Declination (J2000)')
    fig1.axis_labels.hide()
    # No tickmarks
    fig1.tick_labels.hide()
    # ticks colors
    fig0.ticks.set_color('black') 
    fig1.ticks.set_color('black') 

    fig0.add_label(0.95, 0.95, r"N$_2$H$^+$(1-0)", relative=True,
                   horizontalalignment='right')
    fig1.add_label(0.95, 0.95, r"NH$_3$(1,1)", relative=True,
                   horizontalalignment='right')

    fig0.add_colorbar(box=bar_pos0, box_orientation='horizontal',
                      ticks=[0.05, 0.15, 0.25],
                      axis_label_text=r"(km s$^{-1}$)")
    fig1.add_colorbar(box=bar_pos1, box_orientation='horizontal',
                      ticks=[0.05, 0.15, 0.25],
                      axis_label_text=r"(km s$^{-1}$)")
    fig0.set_system_latex(True)
    fig1.set_system_latex(True)
    add_markers_source(fig0, yso_color='yellow', cond_color='0.7')
    add_markers_source(fig1, yso_color='yellow', cond_color='0.7')
    #
    fig.savefig('figures/B5_sigma_v_maps.pdf', bbox_inches='tight')

if do_compare_V:
    plt.ion()
    plt.close()
    V_NH3 = fits.getdata(mergedFile_NH3_Vlsr)
    eV_NH3 = fits.getdata(mergedFile_NH3_eVlsr)
    V_N2Hp = fits.getdata(mergedFile_N2Hp_Vlsr)
    eV_N2Hp = fits.getdata(mergedFile_N2Hp_eVlsr)
    gd = np.isfinite(V_NH3 * V_N2Hp)
    xrange = np.array([9.93, 10.52])
    #
    xx, yy, KDE_vlsr = my_KDE(V_NH3[gd], V_N2Hp[gd],
                              weights=1./(eV_NH3[gd]*eV_N2Hp[gd]),
                              xmin=xrange[0], xmax=xrange[1],
                              ymin=xrange[0], ymax=xrange[1])
    fig, ax = plt.subplots(figsize=(4, 4))
    levels = levels_norm * KDE_vlsr.max()
    cfset = ax.contourf(xx, yy, KDE_vlsr, colors=color_levels, zorder=2,
                        levels=levels_norm_l * KDE_vlsr.max())
    plt.xlabel(r"$V_{\rm LSR}$(NH$_3$) (km s$^{-1}$)")
    plt.ylabel(r"$V_{\rm LSR}$(N$_2$H$^+$) (km s$^{-1}$)")
    ax.set_xlim(xrange)
    ax.set_ylim(xrange)
    ax.plot(xrange, xrange, color='k', zorder=10)
    ax.plot(xrange, xrange+0.05, color='k', ls=':', zorder=11)
    ax.plot(xrange, xrange-0.05, color='k', ls=':', zorder=12)
    #
    ax.text(0.72, 0.26, r'Barnard 5', horizontalalignment='center',
            transform=ax.transAxes)
    ax.text(0.72, 0.20, r'Ions vs Neutrals', horizontalalignment='center',
            transform=ax.transAxes)
    ax.text(0.72, 0.14, r'Centroid velocity, $V_{\rm LSR}$',
            horizontalalignment='center',
            transform=ax.transAxes)
    #
    ax.text(10.35, 10.42, r'$\Delta V_{\rm LSR}$ = $\pm$0.05 km s$^{-1}$',
            horizontalalignment='right')
    ax.annotate(text='', xy=(10.40, 10.35), xytext=(10.40, 10.45),
                arrowprops=dict(arrowstyle='<->'))
    ax.set_xticks([10.0, 10.2, 10.4])
    ax.set_yticks([10.0, 10.2, 10.4])
    fig.savefig('figures/B5_compare_Vlsr_KDE.pdf', bbox_inches='tight')

if do_compare_dv:
    plt.ion()
    plt.close()
    dv_NH3 = fits.getdata(mergedFile_NH3_sigma)
    edv_NH3 = fits.getdata(mergedFile_NH3_esigma)
    dv_N2Hp = fits.getdata(mergedFile_N2Hp_sigma)
    edv_N2Hp = fits.getdata(mergedFile_N2Hp_esigma)
    gd = np.isfinite(dv_NH3 * dv_N2Hp * edv_N2Hp * edv_NH3)
    xrange = np.array([0.045, 0.235])
    #
    dv_diff = np.round(np.mean(dv_N2Hp[gd] - dv_NH3[gd]), decimals=3)
    dv_diff = 0.015
    print('Mean difference = {0:.3}'.format(dv_diff))
    xx, yy, KDE_dv = my_KDE(dv_NH3[gd], dv_N2Hp[gd],
                            weights=1./(edv_NH3[gd]*edv_N2Hp[gd]),
                            xmin=xrange[0], xmax=xrange[1],
                            ymin=xrange[0], ymax=xrange[1])
    fig, ax = plt.subplots(figsize=(4, 4))
    levels = levels_norm * KDE_dv.max()
    cfset = ax.contourf(xx, yy, KDE_dv, colors=color_levels, zorder=2,
                        levels=levels_norm_l * KDE_dv.max())
    plt.xlabel(r"$\sigma_{v}$(NH$_3$) (km s$^{-1}$)")
    plt.ylabel(r"$\sigma_{v}$(N$_2$H$^+$) (km s$^{-1}$)")
    ax.set_xlim(xrange)
    ax.set_ylim(xrange)
    ax.plot(xrange, xrange, color='k', zorder=10)
    ax.plot(xrange, xrange + dv_diff, color='k', zorder=11, ls=':')
    ax.set_xticks([0.05, 0.10, 0.15, 0.2])
    ax.set_yticks([0.05, 0.10, 0.15, 0.2])
    ax.text(0.72, 0.26, r'Barnard 5', horizontalalignment='center',
            transform=ax.transAxes)
    ax.text(0.72, 0.20, r'Ions vs Neutrals', horizontalalignment='center',
            transform=ax.transAxes)
    ax.text(0.72, 0.14, r'Velocity dispersion, $\sigma_v$',
            horizontalalignment='center',
            transform=ax.transAxes)
    ax.text(0.72, 0.08, r'Thermal + Non-Thermal',
            horizontalalignment='center',
            transform=ax.transAxes)
    #
    label_text = r'$\Delta \sigma_v$ = '
    label_text += '{0:.3}'.format(dv_diff)
    label_text += ' km s$^{-1}$'
    ax.text(0.195, 0.2 + dv_diff, label_text, horizontalalignment='right')
    ax.annotate(text='', xy=(0.2, 0.2), xytext=(0.2, 0.2 + dv_diff),
                arrowprops=dict(arrowstyle='<->'))
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
    fig = plt.figure(figsize=(4, 3.7))
    fig0 = aplpy.FITSFigure(fits.PrimaryHDU(ratio_dv, header=hd), figure=fig)
    fig0.show_colorscale(cmap='RdYlBu', vmin=0.5, vmax=1.5)
    fig0.set_nan_color(NaN_color)
    fig0.show_contour(file_NH3_11_match_TdV, levels=clev_NH3_11, colors='k',
                      linewidths=0.5)
    fig0.axis_labels.set_xtext('Right Ascension (J2000)')
    fig0.axis_labels.set_ytext('Declination (J2000)')
    fig0.ticks.set_color('black')
    fig0.add_beam(color='black')
    fig0.add_scalebar(scale_bar)
    fig0.scalebar.set_label(scale_text)
    colorbar_text = r"$\sigma_v$(NH$_3$)/$\sigma_v$(N$_2$H$^+$)"
    fig0.add_colorbar(box=[0.50, 0.8, 0.25, 0.05], box_orientation='horizontal',
                      ticks=[0.5, 1, 1.5],
                      axis_label_text=colorbar_text)
    fig0.set_system_latex(True)
    add_markers_source(fig0, yso_color='yellow', cond_color='0.7')
    fig.savefig('figures/B5_ratio_sigma_v_maps.pdf', bbox_inches='tight')

if do_compare_dv_nt:
    plt.ion()
    plt.close()
    tk_NH3 = fits.getdata(mergedFile_NH3_Tk)*u.K
    dv_NH3 = fits.getdata(mergedFile_NH3_sigma)*u.km/u.s
    edv_NH3 = fits.getdata(mergedFile_NH3_esigma)*u.km/u.s
    dv_N2Hp = fits.getdata(mergedFile_N2Hp_sigma)*u.km/u.s
    edv_N2Hp = fits.getdata(mergedFile_N2Hp_esigma)*u.km/u.s
    gd = np.isfinite(dv_NH3 * dv_N2Hp * edv_N2Hp * edv_NH3 / tk_NH3)
    xrange = np.array([0.0, 0.22])
    #
    # try the KDE
    dv_nt_NH3 = np.sqrt(dv_NH3[gd]**2 - sigma_thermal(17., tk=tk_NH3[gd])**2).to(u.km/u.s)
    dv_nt_N2Hp = np.sqrt(dv_N2Hp[gd]**2 - sigma_thermal(29., tk=tk_NH3[gd])**2).to(u.km/u.s)
    cs_gas = sigma_thermal(2.37, tk=tk_NH3[gd]).to(u.km/u.s)
    Mach_NH3 = (dv_nt_NH3 / cs_gas).decompose().value
    Mach_N2Hp = (dv_nt_N2Hp / cs_gas).decompose().value
    #
    wt = 1. / (edv_NH3[gd] * edv_N2Hp[gd])
    gd = np.isfinite(dv_nt_N2Hp * dv_nt_NH3 * wt)
    #
    mean_sigma_nt_ion = np.nanmean(dv_nt_N2Hp[gd])
    mean_sigma_nt_neutral = np.nanmean(dv_nt_NH3[gd])
    mean_sigmasigma = np.nanmean(np.sqrt(dv_nt_NH3[gd] * dv_nt_N2Hp[gd]))
    p_range = np.array([0.15865 * 100, 50., 0.84135 * 100])
    #
    sigma_mean_nt_ion = np.round(np.percentile(dv_nt_N2Hp[gd], p_range), decimals=3)
    sigma_mean_nt_neutral = np.round(np.percentile(dv_nt_NH3[gd], p_range), decimals=3)
    label_sigma_nt_NH3 = '<sigma_nt,NH3> = {0:.3f} _{1:.3f} ^+{2:5.3f}'
    print(label_sigma_nt_NH3.format(sigma_mean_nt_neutral[1],
                                    sigma_mean_nt_neutral[0] - sigma_mean_nt_neutral[1],
                                    sigma_mean_nt_neutral[2] - sigma_mean_nt_neutral[1]))
    label_sigma_nt_N2Hp = '<sigma_nt,N2Hp> = {0:.3f} _{1:.3f} ^+{2:5.3f}'
    print(label_sigma_nt_N2Hp.format(sigma_mean_nt_ion[1],
                                     sigma_mean_nt_ion[0] - sigma_mean_nt_ion[1],
                                     sigma_mean_nt_ion[2] - sigma_mean_nt_ion[1]))
    print('<sqrt(sigma_nt,NH3 * sigma_nt,N2H+)> = {0:.3}'.format(mean_sigmasigma))
    #
    xx, yy, KDE_dv = my_KDE(dv_nt_NH3[gd].value, dv_nt_N2Hp[gd].value,
                            weights=wt[gd].value,
                            xmin=xrange[0], xmax=xrange[1],
                            ymin=xrange[0], ymax=xrange[1])
    cs_tk_mean = sigma_thermal(2.37, tk=9.7*u.K).to(u.km/u.s).value
    fig, ax = plt.subplots(figsize=(4, 4))
    levels = levels_norm * KDE_dv.max()
    cfset = ax.contourf(xx, yy, KDE_dv, colors=color_levels, zorder=2,
                        levels=levels_norm_l * KDE_dv.max())
    plt.xlabel(r"$\sigma_{NT}$(NH$_3$) (km s$^{-1}$)")
    plt.ylabel(r"$\sigma_{NT}$(N$_2$H$^+$) (km s$^{-1}$)")
    ax.set_xlim(xrange)
    ax.set_ylim(xrange)
    ax.plot(xrange, xrange, color='k', zorder=10)
    ax.set_xticks([0, 0.05, 0.10, 0.15, 0.2])
    ax.set_yticks([0, 0.05, 0.10, 0.15, 0.2])
    ax.text(0.72, 0.26, r'Barnard 5', horizontalalignment='center',
            transform=ax.transAxes)
    ax.text(0.72, 0.20, r'Ions vs Neutrals', horizontalalignment='center',
            transform=ax.transAxes)
    ax.text(0.72, 0.14, r'Velocity dispersion, $\sigma_{NT}$',
            horizontalalignment='center',
            transform=ax.transAxes)
    ax.text(0.72, 0.08, r'Non-Thermal',
            horizontalalignment='center',
            transform=ax.transAxes)
    #
    ax.plot(xrange, xrange + 0.03, color='k', zorder=11, ls=':')
    ax.text(0.175, 0.205, r'$\Delta \sigma_v$ = 0.03 km s$^{-1}$',
            horizontalalignment='right')
    ax.annotate(text='', xy=(0.18, 0.18), xytext=(0.18, 0.21),
                arrowprops=dict(arrowstyle='<->'))
    ax.plot(cs_tk_mean*np.array([0.5, 1]), cs_tk_mean * np.array([1, 1]),
            ls='--', color='#80b1d3')
    ax.plot(cs_tk_mean * np.array([0.5, 0]), cs_tk_mean * np.array([1, 1]) * 0.5,
            ls='--', color='#80b1d3')
    ax.plot(cs_tk_mean * np.array([1, 1]), cs_tk_mean*np.array([0.5, 1]),
            ls='--', color='#80b1d3')
    ax.plot(cs_tk_mean * np.array([1, 1]) * 0.5, cs_tk_mean*np.array([0.5, 0]),
            ls='--', color='#80b1d3')
    ax.text(0.085, 0.025, r'$\mathcal{M}_s=0.5$',
            horizontalalignment='right')
    ax.text(0.085, 0.180, r'$\mathcal{M}_s=1$',
            horizontalalignment='right')
    fig.savefig('figures/B5_compare_sigma_NT_KDE.pdf', bbox_inches='tight')

    fig2, ax2 = plt.subplots(figsize=(4, 4))
    xrange2 = [0., 1.3]
    gd_N2Hp = np.isfinite(Mach_N2Hp)
    gd_NH3 = np.isfinite(Mach_NH3)
    KDE_Mach_N2Hp = stats.gaussian_kde(Mach_N2Hp[gd_N2Hp])
    KDE_Mach_NH3 = stats.gaussian_kde(Mach_NH3[gd_NH3])
    positions = np.arange(xrange2[0], xrange2[1], 0.01)
    ax2.plot(positions, KDE_Mach_N2Hp(positions), color=color_N2Hp, label=r'N$_2$H$^+$')
    ax2.plot(positions, KDE_Mach_NH3(positions), color=color_NH3, label=r'NH$_3$')
    #
    ax2.text(0.78, 0.8,  r'N$_2$H$^+$', color=color_N2Hp, transform=ax2.transAxes, size=14, weight=60)
    ax2.text(0.78, 0.87,  r'NH$_3$', color=color_NH3, transform=ax2.transAxes, size=14, weight=60)
    ax2.set_xlabel(r"Sonic Mach number, $\mathcal{M}_{s}$")
    ax2.set_ylabel(r"Density")
    ax2.set_xlim(xrange2)
    ax2.set_xticks([0, 0.3, 0.6, 0.9, 1.2])
    Mach_N2Hp_med = np.round(np.median(Mach_N2Hp[gd_N2Hp]), decimals=2)
    Mach_NH3_med = np.round(np.median(Mach_NH3[gd_NH3]), decimals=2)
    #
    sigma_Mach_N2Hp = np.round(np.percentile(Mach_N2Hp[gd_N2Hp], p_range), decimals=2)
    sigma_Mach_NH3 = np.round(np.percentile(Mach_NH3[gd_NH3], p_range), decimals=2)
    label_Mach_NH3 = '<Mach,NH3> = {0:.2f} _{1:.2f} ^+{2:.2f}'
    label_Mach_N2Hp = '<Mach,N2Hp> = {0:.2f} _{1:.2f} ^+{2:.2f}'
    print(label_Mach_N2Hp.format(sigma_Mach_N2Hp[1],
                                 sigma_Mach_N2Hp[0] - sigma_Mach_N2Hp[1],
                                 sigma_Mach_N2Hp[2] - sigma_Mach_N2Hp[1]))
    print(label_Mach_NH3.format(sigma_Mach_NH3[1],
                                sigma_Mach_NH3[0] - sigma_Mach_NH3[1],
                                sigma_Mach_NH3[2] - sigma_Mach_NH3[1]))
    ax2.plot([Mach_N2Hp_med, Mach_N2Hp_med], [0, 0.25], color=color_N2Hp, ls='--')
    ax2.plot([Mach_NH3_med, Mach_NH3_med], [0, 0.5], color=color_NH3, ls='--')
    ax2.text(Mach_N2Hp_med-0.07, 0.28, r'{0}$_{{{1:.2f}}}^{{+{2:.2f}}}$'.format(Mach_N2Hp_med,
                                                        sigma_Mach_N2Hp[0] - sigma_Mach_N2Hp[1],
                                                        sigma_Mach_N2Hp[2] - sigma_Mach_N2Hp[1]),
             horizontalalignment='left', color=color_N2Hp,
             size=13, weight=60)
    #r'{0}'.format(Mach_NH3_med),
    ax2.text(Mach_NH3_med-0.07, 0.51, r'{0}$_{{{1:.2f}}}^{{+{2:.2f}}}$'.format(Mach_NH3_med,
                                                        sigma_Mach_NH3[0] - sigma_Mach_NH3[1],
                                                        sigma_Mach_NH3[2] - sigma_Mach_NH3[1]),
             horizontalalignment='left', color=color_NH3,
             size=13, weight=60)
    fig2.savefig('figures/B5_compare_sigma_NT_KDE_dist.pdf', bbox_inches='tight')

