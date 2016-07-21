# Licensed under an MIT open source license - see LICENSE


import numpy as np
import sys
from datetime import datetime
from itertools import repeat


from turbustat.statistics import stats_wrapper, statistics_list


from analysis_funcs import (load_and_reduce, single_input, sort_distances,
                            files_sorter)

np.random.seed(248954785)


def run_all(fiducial, simulation_runs, statistics, savename,
            pool=None, verbose=True,
            multi_timesteps=False, noise_added=False):
    '''
    Given a fiducial set and a series of sets to compare to, loop
    through and compare all sets and their time steps. Return an array of
    the distances.

    Parameters
    ----------
    verbose : bool, optional
        Prints out the time when completing a set.
    multi_timesteps : bool, optional
        If multiple timesteps are given for each simulation run, parallelize
        over the timesteps. If only one is given, parallelize over the
        simulation runs.
    '''

    if verbose:
        # print "Simulation runs to be analyzed: %s" % (simulation_runs)
        print "Started at " + str(datetime.now())

    if multi_timesteps:
        # Distances will be stored in an array of dimensions
        # # statistics x # sim runs x # timesteps
        distances_storage = np.zeros((len(statistics),
                                      len(simulation_runs),
                                      len(fiducial)))

        print distances_storage.shape

        for i, key in enumerate(simulation_runs.keys()):
            timesteps = simulation_runs[key]

            if verbose:
                print "On Simulation %s/%s" % (i + 1, len(simulation_runs))
                print str(datetime.now())
            if pool is not None:

                distances = pool.map(single_input, zip(fiducial,
                                                       timesteps,
                                                       repeat(statistics),
                                                       repeat(noise_added)))

                # If there aren't the maximum number of timesteps, pad the
                # output to match the max.
                if len(distances) < len(fiducial):
                    diff = len(fiducial) - len(distances)
                    for d in range(diff):
                        distances.append(dict.fromkeys(statistics, np.NaN))

                distances_storage[:, i, :] = \
                    sort_distances(statistics, distances).T

            else:
                for ii, timestep in enumerate(timesteps):
                    fiducial_dataset = load_and_reduce(fiducial[ii])
                    testing_dataset = load_and_reduce(timestep)
                    if i == 0:
                        distances, fiducial_models = \
                            stats_wrapper(fiducial_dataset, testing_dataset,
                                          statistics=statistics)
                        all_fiducial_models = fiducial_models
                    else:
                        distances = \
                            stats_wrapper(fiducial_dataset, testing_dataset,
                                          fiducial_models=all_fiducial_models,
                                          statistics=statistics)
                    distances = [distances]
                    distances_storage[:, i, ii:ii + 1] = \
                        sort_distances(statistics, distances).T

    else:
        distances_storage = np.zeros((len(statistics),
                                      len(simulation_runs)))

        if pool is not None:
            # print zip(repeat(fiducial), simulation_runs.values(),
            #           repeat(statistics))
            # print blah
            distances = pool.map(single_input, zip(repeat(fiducial),
                                                   simulation_runs.values(),
                                                   repeat(statistics),
                                                   repeat(noise_added)))

            distances_storage = sort_distances(statistics, distances).T

        else:
            # Load the fiducial in
            fiducial_dataset = load_and_reduce(fiducial)

            for i, key in enumerate(simulation_runs.keys()):
                testing_dataset = \
                    load_and_reduce(simulation_runs[key])
                if i == 0:
                    distances, fiducial_models = \
                        stats_wrapper(fiducial_dataset, testing_dataset,
                                      statistics=statistics)
                    all_fiducial_models = fiducial_models
                else:
                    distances = \
                        stats_wrapper(fiducial_dataset, testing_dataset,
                                      fiducial_models=all_fiducial_models,
                                      statistics=statistics)
                distances = [distances]
                distances_storage[:, i:i + 1] = \
                    sort_distances(statistics, distances).T

    return distances_storage

if __name__ == "__main__":

    # Call as:
    # python output.py path/to/folder/ 0 0 1 max fiducial0 T
    #  /lustre/home/ekoch/results/
    # The args correspond to: directory, fiducial number, face,
    # comparison face, time steps to use, output file prefix,
    # use multiple cores?, save_direc

    # statistics = ["Wavelet", "MVC", "PSpec", "Bispectrum", "DeltaVariance",
    #               "Genus", "VCS", "VCA", "Tsallis", "PCA", "SCF", "Cramer",
    #               "Skewness", "Kurtosis", "VCS_Density", "VCS_Velocity",
    #               "PDF_Hellinger", "PDF_KS", "Dendrogram_Hist",
    #               "Dendrogram_Num"]

    statistics = statistics_list

    statistics = ["Wavelet"] #, "MVC", "PSpec", "Bispectrum", "SCF", "VCA"]

    print "Statistics to run: %s" % (statistics)
    num_statistics = len(statistics)

    # Read in cmd line args

    # Read in all files in the given directory

    print("Running with parameters:")
    print("Data path: %s" % sys.argv[1])
    print("Fiducial %s" % sys.argv[2])
    print("Face 1 %s" % sys.argv[3])
    print("Face 2 %s" % sys.argv[4])
    print("Timesteps %s" % sys.argv[5])
    print("Save name %s" % sys.argv[6])
    print("Run in parallel? %s" % sys.argv[7])
    print("Add noise? %s" % sys.argv[8])
    print("Output to: %s" % sys.argv[9])

    PREFIX = str(sys.argv[1])

    try:
        fiducial_num = int(sys.argv[2])
    except ValueError:
        fiducial_num = str(sys.argv[2])
    face = int(sys.argv[3])
    comp_face = int(sys.argv[4])
    try:
        timesteps = int(sys.argv[5])
    except ValueError:
        timesteps = str(sys.argv[5])
    save_name = str(sys.argv[6])
    MULTICORE = str(sys.argv[7])
    if MULTICORE == "T":
        MULTICORE = True
    else:
        MULTICORE = False
    noise_added = str(sys.argv[8])
    if noise_added == "T":
        noise_added = True
    else:
        noise_added = False
    output_direc = str(sys.argv[9])

    # Run on single timesteps when each simulation is at ~ one free fall time
    if timesteps == "freefall":
        # These are the output times that correspond to ~ free-fall time
        des_tsteps = \
            np.array([25, 26, 21, 21, 23, 25, 21, 21, 25, 26, 21, 22, 23, 26,
                      21, 21, 30, 30, 27, 28, 30, 30, 27, 23, 30, 30, 27, 25,
                      30, 30, 24, 25])

        timesteps = dict.fromkeys(np.arange(32))
        for key, val in zip(timesteps.keys(), des_tsteps):
            timesteps[key] = val

        # Make sure this matches the number of design sims
        assert len(timesteps) == 32

        fid_tsteps = {0: 25, 1: 24, 2: 26, 3: 26, 4: 26}
        assert len(fid_tsteps) == 5
    else:
        fid_tsteps = None

    # Sigma for COMPLETE NGC1333 data using signal-id (normal dist)
    # Note that the mean is forced to 0
    # rms_noise = 0.1277369117707014 / 2.  # in K

    # Trying noise levels scaled by their brightness distribs
    # Set whether we have multiple timesteps for each set
    condition_multi_tsteps = timesteps is 'last' or \
        isinstance(timesteps, dict)
    if condition_multi_tsteps:
        multi_timesteps = False
    else:
        multi_timesteps = True

    fiducials, designs, timesteps_labels = \
        files_sorter(PREFIX, timesteps=timesteps,
                     append_prefix=True, fiducial_timesteps=fid_tsteps)

    if MULTICORE:

        use_mpi = False
        if use_mpi:
            from mpipool import MPIPool
            pool = MPIPool(loadbalance=False)

            if not pool.is_master():
                # Wait for instructions from the master process.
                pool.wait()
                sys.exit(0)
        else:
            from multiprocessing import Pool
            # Default to 10 for now. Will change if this works.
            # pool = Pool(processes=12)
            pool = Pool(processes=3)
    else:
        pool = None

    if fiducial_num == "fid_comp":  # Run all the comparisons of fiducials

        print "Fiducials to compare %s" % (fiducials[face].keys())
        fiducial_index = []
        fiducial_col = []

        # number of comparisons b/w all fiducials
        num_comp = (len(fiducials[face])**2. - len(fiducials[face])) / 2
        # Change dim 2 to match number of time steps
        if isinstance(fiducials[0][0], list):
            distances_storage = np.zeros((num_statistics, num_comp,
                                          len(fiducials[0][0])))
        else:
            distances_storage = np.zeros((num_statistics, num_comp, 1))
        posn = 0
        prev = 0
        # no need to loop over the last one
        for fid_num, i in zip(fiducials[face].keys()[:-1],
                              np.arange(len(fiducials[comp_face]) - 1, 0, -1)):
            posn += i
            comparisons = fiducials[comp_face].copy()

            for key in range(fid_num + 1):
                del comparisons[key]
            partial_distances = \
                run_all(fiducials[face][fid_num], comparisons,
                        statistics, save_name, pool=pool,
                        noise_added=noise_added,
                        multi_timesteps=multi_timesteps, verbose=True)
            shape = distances_storage[:, prev:posn, :].shape
            distances_storage[:, prev:posn, :] = partial_distances.reshape(shape)
            prev += i

            fiducial_index.extend(fiducials[comp_face].keys()[fid_num + 1:])

            fiducial_col.extend(
                [posn - prev] * len(fiducials[comp_face].keys()[fid_num:]))

        # consistent naming with non-fiducial case
        simulation_runs = fiducial_index
        # face = comp_face
    else:  # Normal case of comparing to single fiducial

        distances_storage = \
            run_all(fiducials[face][fiducial_num],
                    designs[comp_face], statistics, save_name,
                    pool=pool,
                    multi_timesteps=multi_timesteps,
                    noise_added=noise_added)

        simulation_runs = designs[comp_face].keys()
        fiducial_index = [fiducial_num] * len(designs.keys())

    # If using timesteps 'max', some comparisons will remain zero
    # To distinguish a bit better, set the non-comparisons to zero
    distances_storage[np.where(distances_storage == 0)] = np.NaN

    filename = save_name + "_fiducial" + str(fiducial_num) + "_" + \
        str(face) + "_" + str(comp_face) + "_distance_results.h5"

    from pandas import DataFrame, HDFStore, concat

    # Save data for each statistic in a dataframe.
    # Each dataframe is saved in a single hdf5 file

    store = HDFStore(output_direc + filename)

    for i in range(num_statistics):
        # If timesteps is 'max', there will be different number of labels
        # in this case, don't bother specifying column names.
        # This also applies for when timesteps is given as free-fall.
        if 'max' in timesteps or isinstance(timesteps, dict):
            df = DataFrame(distances_storage[i], index=simulation_runs)
        else:
            df = DataFrame(distances_storage[i], index=simulation_runs,
                           columns=timesteps_labels[0][face])

        # if not "Fiducial" in df.columns:
        #    df["Fiducial"] = Series(fiducial_index, index=df.index)
        if statistics[i] in store:
            existing_df = store[statistics[i]]
            if len(existing_df.index) == len(df.index):
                store[statistics[i]] = df
            else:  # Append on
                for ind in df.index:
                    if ind in list(existing_df.index):
                        existing_df.ix[ind] = df.ix[ind]
                    else:
                        existing_df = concat([existing_df, df])
                    store[statistics[i]] = existing_df
        else:
            store[statistics[i]] = df

    store.close()

    if MULTICORE:
        pool.close()

    print "Done at " + str(datetime.now())
