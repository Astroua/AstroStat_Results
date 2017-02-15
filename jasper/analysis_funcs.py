
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

    if noise_added:
        vcs_break = -0.5
    else:
        vcs_break = -0.8

    # Find the minimum intensity values for computing dendrograms
    fid_noise = 0.1 * np.nanpercentile(fiducial_dataset["cube"][0], 98)
    test_noise = 0.1 * np.nanpercentile(testing_dataset["cube"][0], 98)

    dendro_params_fid = {"min_value": 2 * fid_noise, "min_npix": 10}
    dendro_params_test = {"min_value": 2 * test_noise, "min_npix": 10}
    dendro_params = [dendro_params_fid, dendro_params_test]

    # The simulations (in set 8) are not benefiting from adding a break,
    # and including one makes the fitting for a few unstable
    vca_break = None

    distances = stats_wrapper(fiducial_dataset, testing_dataset,
                              statistics=statistics, multicore=True,
                              vca_break=vca_break, vcs_break=vcs_break,
                              dendro_params=dendro_params)
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
                 fiducial_identifier="Fiducial"):
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

        if not identified:
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
