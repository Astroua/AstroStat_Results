
import numpy as np
from astropy.io.fits import getdata
from pandas import DataFrame, read_csv
from pandas import HDFStore
from itertools import combinations, repeat
from datetime import datetime
import warnings

from wrapping_function import stats_wrapper
from analysis_funcs import files_sorter
from turbustat.statistics import statistics_list

'''
COMPLETE comparisons.

Inter-compare the observational cubes, as well as running against
the Fiducial runs
'''


def obs_to_obs(file_list, statistics, pool=None):
    '''
    Pass a list of the observational cubes.
    Reduce the data, run the statistics and output
    a csv table of the comparisons.

    If a pool is passed, it runs in parallel.
    '''

    num_comp = len(file_list) * (len(file_list) - 1) / 2

    distances = \
        DataFrame([(i.split("/")[-1], j.split("/")[-1]) for i, j in
                   combinations(file_list, 2)],
                  columns=['Fiducial1', 'Fiducial2'])

    dendro_saves = \
        [(i[:-5] + "_dendrostat.pkl",
          j[:-5] + "_dendrostat.pkl")
         for i, j in combinations(file_list, 2)]

    for stat in statistics:
        distances[stat] = np.zeros((num_comp, ))

    generator = zip(combinations(file_list, 2),
                    repeat(statistics),
                    repeat(False),
                    dendro_saves)

    if pool is None:

        for i, combo in enumerate(generator):

            distance_dict = run_comparison(*combo)[0]

            for key in distance_dict.keys():
                distances[key][i] = distance_dict[key]

    else:

        outputs = pool.map(single_input, generator)

        for i, output in enumerate(outputs):

            distance_dict = output[0]

            for key in distance_dict.keys():
                distances[key][i] = distance_dict[key]

    return distances


def obs_to_fid(obs_list, fiducial_dict, statistics, pool=None):
    '''
    Treat observations as the designs.
    '''

    distances = {}

    for stat in statistics:
        distances[stat] = \
            np.zeros((len(obs_list),
                      len(fiducial_dict.keys())))

    for posn, obs in enumerate(obs_list):

        # Give dendrogram save file.
        dendro_saves = [None, obs[:-5] + "_dendrostat.pkl"]
        scf_saves = [None, obs[:-5] + ".scf.pkl"]

        # Create generator
        gen = zip(zip(fiducial_dict.values(), repeat(obs)),
                  repeat(statistics),
                  repeat(True),
                  repeat(dendro_saves),
                  repeat(scf_saves))

        print("On " + str(posn + 1) + "/" + str(len(obs_list)) + " at " +
              str(datetime.now()))

        if pool is not None:
            outputs = pool.map(single_input, gen)
        else:
            outputs = map(single_input, gen)

        for output in outputs:

            pos1 = obs_list.index(output[2])
            pos2 = fiducial_dict.values().index(output[1])

            distance_dict = output[0]

            # Loop through statistics
            for key in distance_dict.keys():
                distances[key][pos1, pos2] = distance_dict[key]

    return distances


def des_to_obs(obs_list, design_dict, statistics, pool=None):
    '''
    Treat observations as the fiducials.
    '''

    distances = {}

    for stat in statistics:
        distances[stat] = \
            np.zeros((len(design_dict.keys()), len(obs_list)))

    for posn, obs in enumerate(obs_list):

        # Give dendrogram save file.
        dendro_saves = [obs[:-5] + "_dendrostat.pkl", None]
        scf_saves = [obs[:-5] + ".scf.pkl", None]

        # Create generator
        gen = zip(zip(repeat(obs), design_dict.values()),
                  repeat(statistics),
                  repeat(True),
                  repeat(dendro_saves),
                  repeat(scf_saves))

        print("On " + str(posn + 1) + "/" + str(len(obs_list)) + " at " +
              str(datetime.now()))

        if pool is not None:
            outputs = pool.map(single_input, gen)
        else:
            outputs = map(single_input, gen)

        for output in outputs:

            pos1 = design_dict.values().index(output[2])
            pos2 = obs_list.index(output[1])

            distance_dict = output[0]

            # Loop through statistics
            for key in distance_dict.keys():
                distances[key][pos1, pos2] = distance_dict[key]

    return distances


def run_comparison(fits, statistics, add_noise, dendro_saves=[None, None],
                   scf_saves=[None, None]):

    fits1, fits2 = fits

    # Derive the property arrays assuming uniform noise (for sims)
    fiducial_dataset = load_and_reduce(fits1)
    testing_dataset = load_and_reduce(fits2)

    noise_ngc1333 = 0.128
    noise_ophA = 0.252
    noise_ic348 = 0.143

    if "SimSuite" in fits1:
        low_freq1 = 4.5 / 128.
        high_freq1 = 16. / 128.

        low_spat1 = 11.
        high_spat1 = 27.
        fid_noise = 0.1 * np.nanpercentile(fiducial_dataset["cube"][0], 98)

        vca_break1 = None

    else:
        if 'ic348' in fits1:
            fid_noise = noise_ic348
            low_spat1 = 19.
            high_spat1 = 53.
        elif 'ngc1333' in fits1:
            fid_noise = noise_ngc1333
            low_spat1 = 14.
            high_spat1 = 26.
        elif "ophA" in fits1:
            fid_noise = noise_ophA
            low_spat1 = 11.5
            high_spat1 = 38.
        else:
            raise ValueError("fits1 is not an expected observational name")

        low_freq1 = 4. / 256.
        high_freq1 = 32. / 256.

        # Frequency cut-offs exclude noise dominate scales
        vca_break1 = None  # -0.7

    if "SimSuite" in fits2:
        low_freq2 = 4.5 / 128.
        high_freq2 = 16. / 128.

        low_spat2 = 11.
        high_spat2 = 27.
        test_noise = 0.1 * np.nanpercentile(testing_dataset["cube"][0], 98)

        vca_break2 = None

    else:
        if 'ic348' in fits2:
            test_noise = noise_ic348
            low_spat2 = 19.
            high_spat2 = 53.
        elif 'ngc1333' in fits2:
            test_noise = noise_ngc1333
            low_spat2 = 14.
            high_spat2 = 26.
        elif "ophA" in fits2:
            test_noise = noise_ophA
            low_spat2 = 11.5
            high_spat2 = 38.
        else:
            raise ValueError("fits1 is not an expected observational name")

        low_freq2 = 4. / 256.
        high_freq2 = 32. / 256.

        # Frequency cut-offs exclude noise dominate scales
        vca_break2 = None  # -0.7

    # The simulations (in set 8) are not benefiting from adding a break in VCA,
    # and including one makes the fitting for a few unstable
    vca_break = [vca_break1, vca_break2]
    vcs_break = -0.5

    noise_value = [fid_noise, test_noise]

    dendro_params_fid = {"min_value": 2 * fid_noise, "min_npix": 10}
    dendro_params_test = {"min_value": 2 * test_noise, "min_npix": 10}
    dendro_params = [dendro_params_fid, dendro_params_test]

    inertial_range = [[low_freq1, low_freq2], [high_freq1, high_freq2]]
    spatial_range = [[low_spat1, low_spat2], [high_spat1, high_spat2]]

    # Set the SCF boundary types.
    if "SimSuite" in fits1:
        periodic_bounds = [True, False]
        vcs_regrid = [100, None]
    elif "SimSuite" in fits2:
        periodic_bounds = [False, True]
        vcs_regrid = [None, 100]
    else:
        periodic_bounds = [False, False]
        vcs_regrid = [None, None]

    distances = stats_wrapper(fiducial_dataset, testing_dataset,
                              statistics=statistics, multicore=True,
                              vca_break=vca_break, vcs_break=vcs_break,
                              dendro_saves=dendro_saves,
                              scf_saves=scf_saves,
                              dendro_params=dendro_params,
                              periodic_bounds=periodic_bounds,
                              noise_value=noise_value,
                              inertial_range=inertial_range,
                              spatial_range=spatial_range,
                              vcs_regrid=vcs_regrid)

    return distances, fits1, fits2


def single_input(a):
    return run_comparison(*a)


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
            list(getdata(prefix_direc + moment_folder +
                         sim_name + file_lab + ".fits", 0, header=True))

        # And the errors
        file_dict[dic_lab + "_error"] = \
            list(getdata(prefix_direc + moment_folder +
                         sim_name + file_lab + ".fits", 1, header=True))

    # There is an issue with the angular cell size in all of the sim headers.
    # This really only matters when comparing to real data.
    # CDELT2 and CDELT1 need to be ~1000 times larger (since they were setup
    # to roughly match the angular scale in the COMPLETE data.)
    if "homeeros" in sim_name:
        for key in file_dict:
            file_dict[key][1]["CDELT1"] *= 1000.
            file_dict[key][1]["CDELT2"] *= 1000.

    return file_dict


def sort_distances(statistics, distances):
    if len(statistics) > 1:
        distance_array = np.empty((len(distances), len(statistics)))
    elif len(statistics) == 1:
        distance_array = np.empty((len(distances), 1))

    for j, dist in enumerate(distances):
        distance_array[j, :] = [dist[stat] for stat in statistics]

    return distance_array


def sort_sim_files(sim_list, sim_labels=np.arange(0, 5),
                   timestep_labels=np.arange(21, 31),
                   sim_type='Fiducial'):
    '''
    Sort by the given labels.
    '''

    sim_dict = dict.fromkeys(sim_labels)

    for key in sim_dict.keys():
        sim_dict[key] = dict.fromkeys(timestep_labels)

    for sim in sim_list:
        for label in sim_labels:
            if sim_type + str(label) + "_" in sim:
                key = label
                break
        else:
            warnings.warn("Cannot find appropriate label for: " + sim)
            continue

        for time in timestep_labels:
            if "_00" + str(time) + "_" in sim:
                tstep = time
                break
        else:
            warnings.warn("Cannot find appropriate timestep for: " + sim)
            continue

        # Remove empty timesteps
        sim_dict[key] =\
            dict((k, v) for k, v in sim_dict[key].iteritems() if v is not None)

        sim_dict[key][tstep] = sim

    return sim_dict


if __name__ == "__main__":

    import os
    import sys

    # statistics =  statistics_list

    # Set to run on the 'good' statistics
    statistics = ["DeltaVariance_Curve", "DeltaVariance_Slope",
                  "DeltaVariance_Centroid_Curve",
                  "DeltaVariance_Centroid_Slope", "VCS_Break",
                  "VCS", "VCS_Large_Scale", "VCS_Small_Scale",
                  "VCA", "PCA", "Cramer", "SCF",
                  "Dendrogram_Hist", "Dendrogram_Num", "MVC"]

    # statistics = ["SCF", "Genus", "DeltaVariance", "Skewness", "Kurtosis"]
    # statistics = ["DeltaVariance"]

    print "Statistics to run: %s" % (statistics)

    obs_dir = sys.argv[1]
    sim_dir = sys.argv[2]
    face = sys.argv[3]

    # Type of comparison
    comparison = str(sys.argv[4])
    valid_comp = ["Obs_to_Fid", "Obs_to_Obs", "Des_to_Obs"]
    if comparison not in valid_comp:
        raise Warning("comparison type give is not valid: " + str(comparison))

    # Set parallel type
    multiproc = sys.argv[5]
    valid_proc = ["MPI", "noMPI", "Single"]
    if multiproc not in valid_proc:
        raise Warning("multiproc type give is not valid: " + str(multiproc))

    # Output results directory
    output_dir = sys.argv[6]

    if output_dir[-1] != "/":
        output_dir += "/"

    save_name = sys.argv[7]

    # Load the list of complete cubes in

    obs_cubes = [obs_dir + f for f in os.listdir(obs_dir) if f[-4:] == 'fits']

    # sim_dir = "/Volumes/RAIDers_of_the_lost_ark/SimSuite8/"

    # Toggle the pool on here

    if multiproc == "MPI":
        from mpipool import MPIPool
        pool = MPIPool(loadbalance=False)

        if not pool.is_master():
            # Wait for instructions from the master process.
            pool.wait()
            sys.exit(0)
    elif multiproc == "noMPI":
        from multiprocessing import Pool
        # pool = Pool(processes=12)
        pool = Pool(processes=4, maxtasksperchild=100)
    else:
        pool = None

    # Do the actual comparisons

    if comparison == "Obs_to_Fid":
        sim_cubes = [sim_dir + f for f in os.listdir(sim_dir)
                     if "Fiducial" in f and "_0" + face + "_" in f]

        sim_dict = sort_sim_files(sim_cubes)

        # Loop through the fiducial sets
        for fid in sim_dict.keys():

            distances = obs_to_fid(obs_cubes, sim_dict[fid], statistics,
                                   pool=pool)

            store = HDFStore(output_dir + save_name + "_fiducial_" + str(fid) +
                             "_face_" + str(face) + ".h5")

            for key in distances.keys():

                df = DataFrame(distances[key],
                               index=[obs.split("/")[-1]
                                      for obs in obs_cubes],
                               columns=[sim.split("/")[-1]
                                        for sim in sim_dict[fid].values()])
                store[key] = df
            store.close()

        if pool is not None:
            pool.close()

    elif comparison == "Des_to_Obs":
        sim_cubes = [sim_dir + f for f in os.listdir(sim_dir)
                     if "Design" in f and "_0" + face + "_" in f]

        sim_dict = sort_sim_files(sim_cubes, sim_labels=np.arange(0, 32),
                                  sim_type="Design")

        # Loop through the fiducial sets
        for des in sim_dict.keys():

            distances = des_to_obs(obs_cubes, sim_dict[des], statistics,
                                   pool=pool)

            store = HDFStore(output_dir + save_name + "_design_" + str(des) +
                             "_face_" + str(face) + ".h5")

            for key in distances.keys():

                df = DataFrame(distances[key],
                               index=[sim.split("/")[-1]
                                      for sim in sim_dict[des].values()],
                               columns=[obs.split("/")[-1]
                                        for obs in obs_cubes])

                store[key] = df

            store.close()

        if pool is not None:
            pool.close()

    else:
        # Pairwise comparisons between the observations only

        complete_distances = obs_to_obs(obs_cubes, statistics, pool=pool)

        if pool is not None:
            pool.close()

        # Now check if there is an existing copy. If so, open it, split out
        # the column names not in the new table, then merge and overwrite.
        filename = os.path.join(output_dir, "complete_comparisons.csv")
        if os.path.exists(filename):
            old_vals = read_csv(filename)

            for col in old_vals.columns:
                if col not in complete_distances.columns:
                    complete_distances[col] = old_vals[col]

        complete_distances.to_csv(filename)

        # for i, stat in enumerate(complete_distances.keys()):

        #     df = DataFrame(complete_distances[stat], index=obs_cubes,
        #                    columns=obs_cubes)

        #     df.to_csv(obs_dir+"complete_comparisons_"+stat+".csv")
