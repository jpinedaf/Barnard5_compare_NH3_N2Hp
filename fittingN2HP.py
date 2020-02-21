#!/usr/bin/env python

import pyspeckit
import astropy.units as u
import astropy.io.fits as fits
import warnings
import numpy as np
import skimage.morphology as immorph
import matplotlib.pyplot as plt
from config import file_N2Hp_base_erode, file_N2Hp_base_erode_rms, \
     thinFile_N2Hp, thickFile_N2Hp, tpeakFile_N2Hp, maskFile_N2Hp

thinFile = thinFile_N2Hp
thickFile = thickFile_N2Hp
tpeakFile = tpeakFile_N2Hp
maskFile = maskFile_N2Hp

freqLine = 93173704000.0 * u.Hz
snr_min = 3.0

ncpu = 4

vmin = 8.8          # mininum velocity of the main hyperfine components
vmax = 12.5         # maximum velocity of the main hyperfine components

xmax = 65.
ymax = 90.

vmin_plot=10.0
vmax_plot=10.5

vmean = 10.2

#
# set some flags

prepareFiles = False

optThin = True
optThick = True

showOptThin = True
showOptThick = True

plot_dv = False

fitSingleSpectrum = False


#
# load the data cube

cube = pyspeckit.Cube(file_N2Hp_base_erode)

#
# these values need to be set

cube.xarr.refX = freqLine
cube.xarr.velocity_convention='radio'
cube.xarr.convert_to_unit('km/s')


if prepareFiles:
    
    rmsMap = fits.getdata(file_N2Hp_base_erode_rms)
    tpeakMap = cube.slice(vmin, vmax, unit='km/s').cube.max(axis=0)

    peakSNR = tpeakMap / rmsMap
    planemask = (peakSNR > snr_min)
    planemask = immorph.opening(planemask, immorph.disk(5))

    hd_cube = cube.header.copy()
    removeKeys = ['NAXIS3', 'CRPIX3', 'CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3', 'SPECSYS']
    for key in removeKeys:
        hd_cube.remove(key)

    hd_cube['WCSAXES'] = 2

    fits.writeto(tpeakFile, tpeakMap, hd_cube, overwrite=True)
    fits.writeto(maskFile, planemask.astype(int), hd_cube, overwrite=True)

else:
    planemask = fits.getdata(maskFile)
    tpeakMap = fits.getdata(tpeakFile)
    rmsMap = fits.getdata(file_N2Hp_base_erode_rms)

    peakSNR = tpeakMap / rmsMap

plt.ion()

cube.Registry.add_fitter('n2hp_vtau', pyspeckit.models.n2hp.n2hp_vtau_fitter, 4)

printOut = []
outFile = []
limMax = []
maxPars = []
fixed = []
verbose = []
guesses = []        # Tex, tau, v_cen, \sigma_v

if optThin:

    printOut.append('Starting optically thin fit...')
    outFile.append(thinFile)
    guesses.append([5.0, 0.1, vmean, 0.2])
    limMax.append([True, False, True, True])
    maxPars.append([250.0, 0.0, vmax, 1.0])
    fixed.append([False, True, False, False])
    verbose.append(1)

if optThick:

    printOut.append('Starting optically thick fit...')   
    outFile.append(thickFile)
    guesses.append([7.0, 2.0, vmean, 0.1])
    limMax.append([True, True, True, True])
    maxPars.append([20.0, 50.0, vmax, 1.0])
    fixed.append([False, False, False, False])
    verbose.append(2)
   
for i in range(len(printOut)):
    print("")
    print('%s' % printOut[i])

    cube.fiteach(fittype='n2hp_vtau', guesses=guesses[i], verbose_level=verbose[i],
                 limitedmin=[True, True, True, True], limitedmax=limMax[i],
                 minpars=[2.75, 0.01, vmin, 0.05], maxpars=maxPars[i], fixed=fixed[i],
                 use_neighbor_as_guess=True, position_order=1/peakSNR, errmap=rmsMap,
                 maskmap=planemask, multicore=ncpu, start_from_point=(46, 58))

    cube.write_fit(outFile[i], overwrite=True)

# 
# some visualization 

cube.Registry.add_fitter('n2hp_vtau', pyspeckit.models.n2hp.n2hp_vtau_fitter, 4)        # is that really needed again??

vmin = []
vmax = []
modelFile = []
plotNames = []

if showOptThin:
    vmin.append(0.05)
    vmax.append(0.25)  
    modelFile.append(thinFile)
    plotNames.append('figures/B5_N2Hp_pyspeckit_poor_thin_fit.png')

if showOptThick:
    vmin.append(0.05)
    vmax.append(0.21)
    modelFile.append(thickFile)
    plotNames.append('figures/B5_N2Hp_pyspeckit_poor_thick_fit.png')


for p, plotName in enumerate(plotNames):

    cube.load_model_fit(modelFile[p], npars=4, npeaks=1, _temp_fit_loc=(xmax,ymax))
    cube.mapplot()
    cube.plot_spectrum(xmax, ymax, plot_fit=True)
    
    if plot_dv:
        cube.mapplot.plane = cube.parcube[3,:,:]
        cube.mapplot(estimator=None, vmin=0.05, vmax=vmax[p])
    else:
        cube.mapplot.plane = cube.parcube[2,:,:]
        cube.mapplot(estimator=None, vmin=vmin_plot, vmax=vmax_plot)

    plt.draw()
    plt.show()
    plt.savefig(plotName)
