
from astropy.io.fits import getdata
import os
import numpy as np
import copy

from wrapping_function import stats_wrapper


'''
Common functions for the different analysis scripts
'''


def timestep_wrapper(fiducial_timestep, testing_timestep, statistics,
                     noise_added):

    # Load the property arrays assuming uniform noise (for sims)
    fiducial_dataset = load_and_reduce(fiducial_timestep)
    testing_dataset = load_and_reduce(testing_timestep)

    # Find the minimum intensity values for computing dendrograms
    fid_noise = 0.1 * np.nanpercentile(fiducial_dataset["cube"][0], 98)
    test_noise = 0.1 * np.nanpercentile(testing_dataset["cube"][0], 98)

    # NOTE: These fit ranges are only valid for the 128^2 simulated cubes!

    # Only fit to the inertial range for power spectrum based statistics
    # Spat. PSpec., VCA, MVC
    inertial_range = [[4.5 / 128.] * 2, [16. / 128.] * 2]

    # Range to fit to delta-variance and wavelets
    # 1/10 to 1/5 of the box avoids the obvious kernel effects on large-scales
    spatial_range = [[11.] * 2, [27.] * 2]

    if noise_added:
        vcs_break = -0.5
        noise_value = [fid_noise, test_noise]
    else:
        vcs_break = -0.8
        noise_value = [-np.inf, -np.inf]

    dendro_params_fid = {"min_value": 2 * fid_noise, "min_npix": 10}
    dendro_params_test = {"min_value": 2 * test_noise, "min_npix": 10}
    dendro_params = [dendro_params_fid, dendro_params_test]

    # We're only comparing sims here, so all spatial bounds are periodic.
    periodic_bounds = [True, True]

    # The simulations (in set 8) are not benefiting from adding a break,
    # and including one makes the fitting for a few unstable
    vca_break = None

    # Spectrally downsample the cubes for VCS (no info on smallest scales)
    vcs_regrid = [True, True]

    distances = stats_wrapper(fiducial_dataset, testing_dataset,
                              statistics=statistics, multicore=True,
                              vca_break=vca_break, vcs_break=vcs_break,
                              vcs_regrid=vcs_regrid,
                              dendro_params=dendro_params,
                              noise_value=noise_value,
                              inertial_range=inertial_range,
                              spatial_range=spatial_range,
                              periodic_bounds=periodic_bounds)
    return distances


def single_input(a):
    return timestep_wrapper(*a)


def load_and_reduce(filename, moment_folder="moments/"):
    '''
    Load the cube in and derive the property arrays.
    '''

    file_dict = {}

    file_labels = ["_moment0", "_centroid", "_linewidth", "_intint"]
    labels = ["moment0", "centroid", "linewidth", "integrated_intensity"]

    # load the cube in
    file_dict['cube'] = list(getdata(filename, header=True))

    prefix_direc = "/".join(filename.split("/")[:-1])
    if len(prefix_direc) != 0:
        prefix_direc = prefix_direc + "/"
    sim_name = os.path.splitext(os.path.basename(filename))[0]

    for dic_lab, file_lab in zip(labels, file_labels):
        file_dict[dic_lab] = \
            list(getdata(os.path.join(prefix_direc, moment_folder,
                                      sim_name + file_lab + ".fits"),
                         0, header=True))

        # And the errors
        file_dict[dic_lab + "_error"] = \
            list(getdata(os.path.join(prefix_direc, moment_folder,
                                      sim_name + file_lab + ".fits"),
                         1, header=True))

    return file_dict


def sort_distances(statistics, distances):
    if len(statistics) > 1:
        distance_array = np.empty((len(distances), len(statistics)))
    elif len(statistics) == 1:
        distance_array = np.empty((len(distances), 1))

    for j, dist in enumerate(distances):
        distance_array[j, :] = [dist[stat] for stat in statistics]

    return distance_array


def files_sorter(folder, fiducial_labels=np.arange(0, 5, 1),
                 design_labels=np.arange(0, 32, 1), timesteps='last',
                 faces=[0, 1, 2], suffix="fits", append_prefix=False,
                 fiducial_timesteps=None, design_identifier="Design",
                 fiducial_identifier="Fiducial", verbose=True):
    '''
    If the entire simulation suite is in one directory, this function
    will spit out appropriate groupings.

    Parameters
    ----------
    folder : str
        Folder where data is.
    fiducial_labels : list or numpy.ndarray, optional
        List of the fiducial numbers.
    design_labels : list or numpy.ndarray, optional
        List of the design numbers.
    timesteps : 'last' or list or numpy.ndarray, optional
        List of timesteps to analyze. If 'last', the last timestep
        found for each simulation is used.
    faces : list
        Faces of the simulations to use.
    suffix : str, optional
        File suffix.
    append_prefix : bool, optional
        Append the complete path to each file.
    fiducial_timesteps : {None, list}, optional
        Optionally provide specific timesteps to use for each fiducial.
    design_identifier : str, optional
        Unique string that identifies a design simulated cube. Default is
        "Design".
    fiducial_identifier : str, optional
        Unique string that identifies a design simulated cube. Default is
        "Fiducial".
    verbose : bool, optional
        Print out some warnings and updates while running.
    '''

    # Get the files and remove any sub-directories.
    files = [f for f in os.listdir(folder) if not os.path.isdir(f) and
             f[-len(suffix):] == suffix]

    # Set up the dictionaries.
    fiducials = dict.fromkeys(faces)
    designs = dict.fromkeys(faces)
    timestep_labels = dict.fromkeys(faces)
    for face in faces:
        fiducials[face] = dict((lab, []) for lab in fiducial_labels)
        designs[face] = dict((lab, []) for lab in design_labels)
        timestep_labels[face] = dict((lab, []) for lab in design_labels)

    # Sort the files
    for f in files:
        # Track if the file was already classified.
        identified = False
        if fiducial_identifier in f:
            for lab in fiducial_labels:
                if not fiducial_identifier + str(lab) + "_" in f:
                    continue
                for face in faces:
                    if "_0" + str(face) + "_" in f:
                        identified = True
                        if append_prefix:
                            fiducials[face][lab].append(os.path.join(folder, f))
                        else:
                            fiducials[face][lab].append(f)

        if design_identifier in f and not identified:
            for lab in design_labels:
                if not design_identifier + str(lab) + "_" in f:
                    continue
                for face in faces:
                    if "_0" + str(face) + "_" in f:
                        identified = True
                        if append_prefix:
                            designs[face][lab].append(os.path.join(folder, f))
                        else:
                            designs[face][lab].append(f)

        if not identified and verbose:
            print("Could not find a category for " + f)

    # Sort and keep only the specified timesteps
    # Can supply different timesteps for the fiducials. This is needed for the
    # common free-fall analysis.
    if fiducial_timesteps is not None:
        _timestep_sort(fiducials, fiducial_timesteps)
    else:
        _timestep_sort(fiducials, timesteps)

    _timestep_sort(designs, timesteps, labels=timestep_labels)

    return fiducials, designs, timestep_labels


def _timestep_sort(d, timesteps, labels=None):
    '''
    Helper function for segmenting by timesteps.
    '''
    for face in d.keys():
        for lab in d[face].keys():
            # Check for empty lists.
            if d[face][lab] == []:
                continue
            d[face][lab].sort()
            if timesteps == 'last':  # Grab the last one
                if labels is not None:
                    labels[face][lab].append(d[face][lab][-1][-16:-14])
                d[face][lab] = d[face][lab][-1]
            elif timesteps == 'max':  # Keep all available
                # Reverse the order so the comparisons are between the highest
                # time steps.
                d[face][lab] = d[face][lab][::-1]
            elif isinstance(timesteps, int):  # Slice out a certain section
                d[face][lab] = d[face][lab][:timesteps]
                if labels is None:
                    continue
                for val in d[face][lab]:
                    labels[face][lab].append(val[-16:-14])
            # Specify a unique timestep for each lab
            elif isinstance(timesteps, dict):
                # if not len(d[face][lab]) == len(timesteps):
                #     print(d[face][lab])
                #     print(timesteps.keys())
                #     raise IndexError("When timesteps is a dict, it must have"
                #                      " the same keys as d.")
                # String to look for in the list of files
                matcher = "_00" + str(timesteps[lab]) + "_"

                match = None
                for f in d[face][lab]:
                    if matcher in f:
                        match = f
                        if labels is not None:
                            labels[face][lab] = timesteps[lab]
                        break
                if match is None:
                    raise Warning("Could not find a match for the timestep"
                                  " {}".format(timesteps[lab]))
                else:
                    d[face][lab] = match
            else:  # Make a copy and loop through the steps
                good_files = copy.copy(d[face][lab])
                for f in d[face][lab]:
                    match = ["_00" + str(step) +
                             "_" in f for step in timesteps]
                    if not any(match):
                        good_files.remove(f)
                    if labels is not None:
                        labels[face][lab].append(f[-16:-14])
                d[face][lab] = good_files


def spectral_regrid_cube(cube, num_chans):
    '''
    Spectrally interpolate a cube to a specified number of channels.
    '''

    new_spec = np.linspace(cube.spectral_extrema[0],
                           cube.spectral_extrema[1], num_chans)

    return cube.spectral_interpolate(new_spec)


def round_up_to_odd(f):
    return np.ceil(f) // 2 * 2 + 1


def make_signal_mask(cube, smooth_chans=200. / 66., min_chan=7, peak_snr=8.,
                     min_snr=5., edge_thresh=0.5):
    '''
    Create a robust signal mask by requiring spatial and spectral
    connectivity.
    '''

    import astropy.units as u
    from astropy.convolution import Box1DKernel
    from signal_id import Noise
    from scipy import ndimage as nd
    from astropy.wcs.utils import proj_plane_pixel_scales
    from astropy.utils.console import ProgressBar
    import skimage.morphology as mo
    import numpy as np
    from radio_beam import Beam
    from itertools import groupby, chain
    from operator import itemgetter
    import matplotlib.pyplot as p

    pixscale = proj_plane_pixel_scales(cube.wcs)[0]

    # # Want to smooth the mask edges
    mask = cube.mask.include()

    # Set smoothing parameters and # consecutive channels.
    smooth_chans = int(round_up_to_odd(smooth_chans))

    # consecutive channels to be real emission.
    num_chans = min_chan

    # Smooth the cube, then create a noise model
    spec_kernel = Box1DKernel(smooth_chans)
    smooth_cube = cube.spectral_smooth(spec_kernel)

    noise = Noise(smooth_cube)
    noise.estimate_noise(spectral_flat=True)
    noise.get_scale_cube()

    snr = noise.snr.copy()

    snr[np.isnan(snr)] = 0.0

    posns = np.where(snr.max(axis=0) >= min_snr)

    bad_pos = np.where(snr.max(axis=0) < min_snr)
    mask[:, bad_pos[0], bad_pos[1]] = False

    # In case single spectra need to be inspected.
    verbose = False

    for i, j in ProgressBar(zip(*posns)):

        # Look for all pixels above min_snr
        good_posns = np.where(snr[:, i, j] > min_snr)[0]

        # Reject if the total is less than connectivity requirement
        if good_posns.size < num_chans:
            mask[:, i, j] = False
            continue

        # Find connected pixels
        sequences = []
        for k, g in groupby(enumerate(good_posns), lambda (i, x): i - x):
            sequences.append(map(itemgetter(1), g))

        # Check length and peak. Require a minimum of 3 pixels above the noise
        # to grow from.
        sequences = [seq for seq in sequences if len(seq) >= 3 and
                     np.nanmax(snr[:, i, j][seq]) >= peak_snr]

        # Continue if no good sequences found
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

        # Final check for the min peak level and ensure all meet the
        # spectral connectivity requirement
        sequences = [seq for seq in sequences if len(seq) >= num_chans and
                     np.nanmax(snr[:, i, j][seq]) >= peak_snr]

        if len(sequences) == 0:
            mask[:, i, j] = False
            continue

        bad_posns = \
            list(set(np.arange(cube.shape[0])) - set(list(chain(*sequences))))

        mask[:, i, j][bad_posns] = False

        if verbose:
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
            p.draw()
            raw_input("Next spectrum?")
            p.clf()

    # initial_mask = mask.copy()

    # Now set the spatial connectivity requirements.

    # The spatial pixel scales in the sim headers are SUPER small
    # choosing this major axis gives an appropriately sized. 5x5 kernel
    beam = Beam(major=1e-3 * u.arcmin)

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

    return masked_cube, noise.scale * masked_cube.unit
