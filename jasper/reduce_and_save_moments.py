
import numpy as np
from spectral_cube import SpectralCube, LazyMask
from astropy import units as u
import os

from turbustat.data_reduction import Mask_and_Moments
from MPI import MPIPool
from multiprocessing import Pool

# Regridding to a common linewidth
from preprocessor import preprocessor
from analysis_funcs import make_signal_mask

'''
Calculate the moments for all of the cubes.
'''


def reduce_and_save(filename, add_noise=False, regrid_linewidth=False,
                    rms_noise=0.001 * u.K, output_path="", cube_output=None,
                    nsig=3, slicewise_noise=True):
    '''
    Load the cube in and derive the property arrays.
    '''

    if add_noise or regrid_linewidth:

        sc = SpectralCube.read(filename)

        if add_noise:
            if rms_noise is None:
                raise TypeError("Must specify value of rms noise.")

            cube = sc.filled_data[:].value

            # Optionally scale noise by 1/10th of the 98th percentile in the
            # cube
            if rms_noise == 'scaled':
                rms_noise = 0.1 * \
                    np.percentile(cube[np.isfinite(cube)], 98) * sc.unit

            from scipy.stats import norm
            if not slicewise_noise:
                cube += norm.rvs(0.0, rms_noise.value, cube.shape)
            else:
                spec_shape = cube.shape[0]
                slice_shape = cube.shape[1:]
                for i in range(spec_shape):
                    cube[i, :, :] += norm.rvs(0.0, rms_noise.value,
                                              slice_shape)

            sc = SpectralCube(data=cube * sc.unit, wcs=sc.wcs,
                              meta={"BUNIT": "K"})

            mask = LazyMask(np.isfinite, sc)
            sc = sc.with_mask(mask)

        if regrid_linewidth:
            # Normalize the cubes to have the same linewidth
            # channels_per_sigma=20 scales to the largest mean line width in
            # SimSuite8 (~800 km/s; Design 22). So effectively everything is
            # "smoothed" to have this line width
            # Intensities are normalized by their 95% value.
            sc = preprocessor(sc, min_intensity=nsig * rms_noise,
                              norm_intensity=True,
                              norm_percentile=95,
                              channels_per_sigma=20)

    else:
        sc = filename

    # Run the same signal masking procedure that was used for the
    # COMPLETE cubes
    if add_noise:
        # The default settings were set based on these cubes
        sc = make_signal_mask(sc)[0]

    reduc = Mask_and_Moments(sc, scale=rms_noise)
    if not add_noise:
        reduc.make_mask(mask=reduc.cube > nsig * reduc.scale)

    reduc.make_moments()
    reduc.make_moment_errors()

    # Remove .fits from filename
    save_name = os.path.splitext(os.path.basename(filename))[0]

    reduc.to_fits(os.path.join(output_path, save_name))

    # Save the noisy cube too
    if add_noise or regrid_linewidth:
        save_name += ".fits"
        if cube_output is None:
            sc.hdu.writeto(os.path.join(output_path, save_name))
        else:
            sc.hdu.writeto(os.path.join(cube_output, save_name))


def single_input(a):
    return reduce_and_save(*a)


if __name__ == "__main__":

    import sys
    import glob
    from itertools import repeat

    folder = str(sys.argv[1])

    output_folder = str(sys.argv[2])

    is_obs = str(sys.argv[3])
    if is_obs == "T":
        is_obs = True
    else:
        is_obs = False

    add_noise = str(sys.argv[4])
    if add_noise == "T":
        add_noise = True
    else:
        add_noise = False

    regrid_linewidth = str(sys.argv[5])
    if regrid_linewidth == "T":
        regrid_linewidth = True
    else:
        regrid_linewidth = False

    if add_noise or regrid_linewidth:
        try:
            cube_output = str(sys.argv[6])
        except IndexError:
            print "Using same output folder for dirty cubes and moments."
            cube_output = output_folder
    else:
        cube_output = output_folder

    # Grab all of the fits files
    fits_files = glob.glob(os.path.join(folder, "*.fits"))

    # Trying noise levels scaled by their brightness distribs
    if add_noise:
        rms_noise = 'scaled'
    elif is_obs:
        rms_noise = None
    else:
        rms_noise = 0.001 * u.K

    np.random.seed(248954785)

    # pool = MPIPool(loadbalance=False)

    # if not pool.is_master():
    #     # Wait for instructions from the master process.
    #     pool.wait()
    #     sys.exit(0)

    pool = Pool(processes=3)

    pool.map(single_input, zip(fits_files,
                               repeat(add_noise),
                               repeat(regrid_linewidth),
                               repeat(rms_noise),
                               repeat(output_folder),
                               repeat(cube_output)))

    pool.close()

    # map(single_input, zip(fits_files,
    #                       repeat(add_noise),
    #                       repeat(regrid_linewidth),
    #                       repeat(rms_noise),
    #                       repeat(output_folder),
    #                       repeat(cube_output)))
