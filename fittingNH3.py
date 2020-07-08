import pyspeckit
import astropy.io.fits as fits
import numpy as np

from spectral_cube import SpectralCube
import astropy.units as u
from skimage.morphology import remove_small_objects, closing, disk, opening

from pyspeckit.spectrum.models import ammonia
from config import file_NH3_11_match, file_NH3_11_match_rms, file_NH3_11_match_TdV,\
                   file_NH3_11_match_mom1, file_NH3_11_match_Tp,\
                   file_NH3_22_match, file_NH3_22_match_rms,\
                   thinFile_NH3, thickFile_NH3, maskFile_NH3

OneOneMom1 = file_NH3_11_match_mom1
OneOnePeak = file_NH3_11_match_Tp
OneOneFile = file_NH3_11_match
RMSFile = file_NH3_11_match_rms
TwoTwoFile = file_NH3_22_match
# Define frequencies used to determine the velocity
snr_min = 8.0

Prepare_files = True

if Prepare_files:
    cube11sc = SpectralCube.read(OneOneFile)
    cube11_v = cube11sc.with_spectral_unit(u.km/u.s,
                        velocity_convention='radio')
    slab = cube11_v.spectral_slab(11.1*u.km/u.s, 9.4*u.km/u.s)

    errmap11 = fits.getdata(RMSFile)
    Tpeak = slab.max(axis=0)
    peaksnr = Tpeak.value/errmap11

    planemask = (peaksnr > snr_min)  # *(errmap11 < 0.15)
    planemask = remove_small_objects(planemask, min_size=15)
    planemask = opening(planemask, disk(5))

    moment1 = slab.moment1( axis=0)

    moment1.write(OneOneMom1, format='fits', overwrite=True)
    Tpeak.write(OneOnePeak, format='fits', overwrite=True)
    hd = fits.getheader(OneOnePeak)
    fits.writeto(maskFile_NH3, planemask.astype(int), hd, overwrite=True)


def b5_cubefit(vmin=9.4, vmax=11.1, do_plot=False, multicore=1,
    file_out='test.fits', tk_mean=10., fix_tk=False):
    """
    Fit NH3(1,1) and (2,2) cubes for B5.
    It fits all pixels with SNR larger than requested. 
    Initial guess is based on moment maps and neighboring pixels. 
    The fitting can be done in parallel mode using several cores, 
    however, this is dangerous for large regions, where using a 
    good initial guess is important. 
    It stores the result in a FITS cube. 

    TODO:
    -convert FITS cube into several FITS files
    -Improve initial guess
    
    Parameters
    ----------
    vmin : numpy.float
        Minimum centroid velocity to plot, in km/s.
    vmax : numpy.float
        Maximum centroid velocity to plot, in km/s.
    do_plot : bool
        If True, then a map of the region to map is shown.
    snr_min : numpy.float
        Minimum signal to noise ratio of the spectrum to be fitted.
    multicore : int
        Numbers of cores to use for parallel processing. 
    """
    errmap11 = fits.getdata(RMSFile)
    Tpeak11 = fits.getdata(OneOnePeak)
    moment1 = fits.getdata(OneOneMom1)
    peaksnr = Tpeak11/errmap11

    planemask = fits.getdata(maskFile_NH3)

    w11 = fits.getdata(file_NH3_11_match_TdV)
    peakloc = np.nanargmax(w11)
    ymax, xmax = np.unravel_index(peakloc, w11.shape)
    
    ## Load FITS files
    cube11 = pyspeckit.Cube(OneOneFile, maskmap=planemask)
    cube22 = pyspeckit.Cube(TwoTwoFile, maskmap=planemask)
    # Stack files
    cubes = pyspeckit.CubeStack([cube11, cube22], maskmap=planemask)
    cubes.unit = "K"
    # Define initial guess
    guesses = np.zeros((6,) + cubes.cube.shape[1:])
    moment1[moment1 < vmin] = vmin+0.2
    moment1[moment1 > vmax] = vmax-0.2
    guesses[0, :, :] = tk_mean  # Kinetic temperature
    guesses[1, :, :] = 7        # Excitation  Temp
    guesses[2, :, :] = 14.5     # log(column)
    guesses[3, :, :] = 0.12     # velocity dispersion
    guesses[4, :, :] = moment1  # Line centroid
    guesses[5, :, :] = 0.5      # F(ortho) - ortho NH3 fraction (fixed)
    if do_plot:
        import matplotlib.pyplot as plt
        plt.imshow(w11*planemask, origin='lower')
        plt.show()
    print('start fit')
    cubes.specfit.Registry.add_fitter('cold_ammonia',
                                      ammonia.cold_ammonia_model(),
                                      6)
    #
    cubes.fiteach(fittype='cold_ammonia',  guesses=guesses,
                  integral=False, verbose_level=3,
                  fixed=[fix_tk, False, False, False, False, True],
                  signal_cut=2,
                  limitedmax=[True, False, False, False, True, True],
                  maxpars=[20, 15, 20, 0.4, vmax, 1],
                  limitedmin=[True, True, True, True, True, True],
                  minpars=[5, 2.8, 12.0, 0.05, vmin, 0],
                  start_from_point=(xmax, ymax),
                  use_neighbor_as_guess=True,
                  position_order=1/peaksnr,
                  errmap=errmap11, multicore=multicore)
    # Store fits into FITS cube
    fitcubefile = fits.PrimaryHDU(
        data=np.concatenate([cubes.parcube, cubes.errcube]),
        header=cubes.header)
    # Prepare FITS header with information about the parameters
    fitcubefile.header.set('PLANE1', 'TKIN')
    fitcubefile.header.set('PLANE2', 'TEX')
    fitcubefile.header.set('PLANE3', 'COLUMN')
    fitcubefile.header.set('PLANE4', 'SIGMA')
    fitcubefile.header.set('PLANE5', 'VELOCITY')
    fitcubefile.header.set('PLANE6', 'FORTHO')
    fitcubefile.header.set('PLANE7', 'eTKIN')
    fitcubefile.header.set('PLANE8', 'eTEX')
    fitcubefile.header.set('PLANE9', 'eCOLUMN')
    fitcubefile.header.set('PLANE10', 'eSIGMA')
    fitcubefile.header.set('PLANE11', 'eVELOCITY')
    fitcubefile.header.set('PLANE12', 'eFORTHO')
    fitcubefile.header.set('CDELT3', 1)
    fitcubefile.header.set('CTYPE3', 'FITPAR')
    fitcubefile.header.set('CRVAL3', 0)
    fitcubefile.header.set('CRPIX3', 1)
    fitcubefile.writeto(file_out, overwrite=True)


#b5_cubefit(multicore=4, file_out=thickFile_NH3, tk_mean=9.5, fix_tk=False)
par_thick = fits.getdata(thickFile_NH3)
Tk = par_thick[0, :, :]
eTk = par_thick[6, :, :]
eTk[eTk == 0] = np.nan
import matplotlib.pyplot as plt
plt.ion()
plt.imshow(Tk, origin='lowest', cmap='inferno', vmin=8, vmax=12.)
plt.colorbar()
plt.contour(eTk, levels=[0.5])

tk_mean = np.round(np.mean(Tk[(eTk < 0.5) * (eTk > 0.)]), decimals=1)
b5_cubefit(multicore=4, file_out=thinFile_NH3, tk_mean=tk_mean, fix_tk=True)
print("Used a Tk = {0} K".format(tk_mean))

par_thick = fits.getdata(thickFile_NH3)
par_thin = fits.getdata(thinFile_NH3)
eTk = par_thick[6, :, :]
mask = (eTk < 0.5) * (eTk > 0.)
par_merged = par_thick*mask + (1 - mask)*par_thin

