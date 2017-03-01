
import astropy.units as u
from astropy.convolution import Box1DKernel
from spectral_cube import SpectralCube
from signal_id import Noise
from scipy import ndimage as nd
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.utils.console import ProgressBar
from astropy.stats import mad_std
import skimage.morphology as mo
import numpy as np
from radio_beam import Beam
from itertools import groupby, chain
from operator import itemgetter
import matplotlib.pyplot as p
from turbustat.data_reduction import Mask_and_Moments
import os

'''
Build a noise model for the COMPLETE cubes and build a good signal mask.
'''


def round_up_to_odd(f):
    return np.ceil(f) // 2 * 2 + 1


def sigma_rob(data, iterations=1, thresh=3.0, axis=None):
    """
    Iterative m.a.d. based sigma with positive outlier rejection.
    """
    noise = mad_std(data, axis=axis)
    for _ in range(iterations):
        ind = (np.abs(data) <= thresh * noise).nonzero()
        noise = mad_std(data[ind], axis=axis)
    return noise


names = ['ic348', 'ngc1333', 'ophA']

for name in names:
    cube = SpectralCube.read('{}.13co.fits'.format(name))

    pixscale = proj_plane_pixel_scales(cube.wcs)[0]

    # cube = cube.with_mask(cube > 3 * scale * u.Jy)

    # # Want to smooth the mask edges
    mask = cube.mask.include()

    # Set smoothing parameters and # consecutive channels.
    smooth_chans = int(round_up_to_odd(200. / 66.))

    # 10 consecutive channels must be above the MAD level to be real emission.
    num_chans = 7

    peak_snr = 4.5
    # Cutoff level
    min_snr = 3
    # Where to cut at the line edges
    edge_thresh = 0.5

    # Smooth the cube, then create a noise model
    spec_kernel = Box1DKernel(smooth_chans)
    smooth_cube = cube.spectral_smooth(spec_kernel)

    noise = Noise(smooth_cube)
    noise.estimate_noise(spectral_flat=True)
    noise.get_scale_cube()

    snr = noise.snr.copy()

    posns = np.where(snr.max(axis=0) >= min_snr)

    bad_pos = np.where(snr.max(axis=0) < min_snr)
    mask[:, bad_pos[0], bad_pos[1]] = False

    # In case single spectra need to be inspected.
    verbose = False

    for i, j in ProgressBar(zip(*posns)):

        spectrum = cube[:, i, j].value

        good_posns = np.where(snr[:, i, j] > min_snr)[0]

        if good_posns.size < num_chans:
            mask[:, i, j] = False
            continue

        sequences = []
        for k, g in groupby(enumerate(good_posns), lambda (i, x): i - x):
            sequences.append(map(itemgetter(1), g))

        # Check length and peak
        # sequences = [seq for seq in sequences if len(seq) >= num_chans and
        #              np.nanmax(smoothed[seq]) >= peak_snr * mad]
        sequences = [seq for seq in sequences if len(seq) >= 5 and
                     np.nanmax(snr[:, i, j][seq]) >= peak_snr]

        if len(sequences) == 0:
            mask[:, i, j] = False
            continue

        # Now take each valid sequence and expand the edges until the smoothed
        # spectrum approaches zero.
        edges = [[seq[0], seq[-1]] for seq in sequences]
        for n, edge in enumerate(edges):
            # Lower side
            if n == 0:
                start_posn = edge[0]
                stop_posn = 0
            else:
                start_posn = edge[0] - edges[n - 1][0]
                stop_posn = edges[n - 1][0]

            for pt in np.arange(start_posn, stop_posn, -1):
                # if smoothed[pt] <= mad * edge_thresh:
                if snr[:, i, j][pt] <= edge_thresh:
                    break

                sequences[n].insert(0, pt)

            # Upper side
            start_posn = edge[1]
            if n == len(edges) - 1:
                stop_posn = cube.shape[0]
            else:
                stop_posn = edges[n + 1][0]

            for pt in np.arange(start_posn, stop_posn, 1):
                # if smoothed[pt] <= mad * edge_thresh:
                if snr[:, i, j][pt] <= edge_thresh:
                    break

                sequences[n].insert(0, pt)

        sequences = [seq for seq in sequences if len(seq) >= num_chans and
                     np.nanmax(snr[:, i, j][seq]) >= peak_snr]

        if len(sequences) == 0:
            mask[:, i, j] = False
            continue

        bad_posns = \
            list(set(np.arange(cube.shape[0])) - set(list(chain(*sequences))))

        mask[:, i, j][bad_posns] = False

        if verbose:
            # p.plot(spectrum)
            # p.plot(smoothed)
            p.subplot(121)
            p.plot(noise.snr[:, i, j])
            p.vlines(np.where(mask[:, i, j])[0][-1], 0,
                     np.nanmax(noise.snr[:, i, j]))
            p.vlines(np.where(mask[:, i, j])[0][0], 0,
                     np.nanmax(noise.snr[:, i, j]))
            p.plot(noise.snr[:, i, j] * mask[:, i, j], 'bD')

            p.subplot(122)
            p.plot(cube[:, i, j], label='Cube')
            p.plot(smooth_cube[:, i, j], label='Smooth Cube')
            p.axvline(np.where(mask[:, i, j])[0][-1])
            p.axvline(np.where(mask[:, i, j])[0][0])
            p.plot(smooth_cube[:, i, j] * mask[:, i, j], 'bD')

            # p.vlines(np.where(mask[:, i, j])[0][-1], 0, np.nanmax(spectrum))
            # p.vlines(np.where(mask[:, i, j])[0][0], 0, np.nanmax(spectrum))
            # p.plot(smoothed * mask[:, i, j], 'bD')
            p.draw()
            raw_input("Next spectrum?")
            p.clf()

    initial_mask = mask.copy()

    # The resolution is about 47 arcsec, but this is just 2 pixels across in
    # the map. I'm going to double this to only include features that are
    # clearly not noise
    beam = Beam(major=2 * 47 * u.arcsec)

    kernel = beam.as_tophat_kernel(pixscale)
    kernel_pix = (kernel.array > 0).sum()

    for i in ProgressBar(mask.shape[0]):
        mask[i] = nd.binary_opening(mask[i], kernel)
        mask[i] = nd.binary_closing(mask[i], kernel)
        mask[i] = mo.remove_small_objects(mask[i], min_size=kernel_pix,
                                          connectivity=2)
        mask[i] = mo.remove_small_holes(mask[i], min_size=kernel_pix,
                                        connectivity=2)

    # Each region must contain a point above the peak_snr
    labels, num = nd.label(mask, np.ones((3, 3, 3)))
    for n in range(1, num + 1):
        pts = np.where(labels == n)
        if np.nanmax(snr[pts]) < peak_snr:
            mask[pts] = False

    masked_cube = cube.with_mask(mask)

    # Save the cube
    masked_cube.write("{}.masked.fits".format(name))

    # Now make
    reduc = Mask_and_Moments(masked_cube, scale=noise.scale)

    reduc.make_moments()
    reduc.make_moment_errors()

    reduc.to_fits(os.path.join("moments/", name))
