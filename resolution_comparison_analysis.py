
from turbustat.data_reduction import Mask_and_Moments
from turbustat.statistics import statistics_list
import os
import sys
import numpy as np
import time
from pandas import DataFrame, read_csv
import statsmodels.api as sm
from astropy.utils.console import ProgressBar
from spectral_cube import SpectralCube
from astropy.io.fits import getheader
from copy import copy
from astropy import units as u

from jasper.wrapping_function import stats_wrapper
from jasper.analysis_funcs import files_sorter, sort_distances

import matplotlib.pyplot as p
import seaborn as sns
sns.set_style("ticks")

p.ioff()

np.random.seed(446784788)

# Note that the moments need to be saved. See reduce_and_save_moments.py
# in jasper/

run_moments = False
# To run on the regridded cubes, enable this
run_on_regridded = True if sys.argv[1] == "T" else False
run_distances = True
run_analysis = True

# path_to_data = "/media/eric/Data_3/Astrostat/SimSuite8/"
path_to_data = os.path.expanduser("~/MyRAID/Astrostat/SimSuite8/")
moments_path = os.path.join(path_to_data, "moments/")

results_path = sys.argv[2]

# Only running on face 0
faces = [0]  # [0, 2]

path_to_256 = os.path.expanduser("~/MyRAID/Astrostat/Fiducial_256/")

if run_on_regridded:
    print("Running on regridded 256 cubes.")

    regrid_path = os.path.join(path_to_256, "regrid")

    # Regrid the 256 cubes to 128 and save them.
    if not os.path.exists(regrid_path):
        os.mkdir(regrid_path)

        fiducials_256, _, _ = \
            files_sorter(path_to_256, append_prefix=True, design_labels=[],
                         fiducial_labels=[256], timesteps=1)

        # Load in a 128 header to regrid to
        hdr = getheader(os.path.join(path_to_data,
                                     "homeerosSimSuite8Fiducial0_flatrho_0021_00_radmc.fits"))

        for face in fiducials_256:
            fid_256 = fiducials_256[face][256][0]

            cube = SpectralCube.read(fid_256)

            regrid_cube = cube.reproject(hdr)

            regrid_cube.write(os.path.join(regrid_path,
                                           fid_256.split("/")[-1]))
    else:
        print("Already running on regridded cubes. Won't regrid again.")

    path_to_256 = os.path.expanduser("~/MyRAID/Astrostat/Fiducial_256/regrid/")
else:
    print("Running on normal 256 cubes.")

moments_256_path = os.path.join(path_to_256, "moments/")

output_path = \
    os.path.join(results_path, "resolution_comparison")
if not os.path.exists(output_path):
    os.mkdir(output_path)


# Sort those from the 128 set, and keep only the first timestep
# B/c the timestep 20 were processed for the fiducials, keep two
# timesteps, then use the last one in the list
fiducials, _, _ = \
    files_sorter(path_to_data, timesteps=2, faces=faces,
                 append_prefix=True, design_labels=[], verbose=False)

# Now the 256 cubes
fiducials_256, _, _ = \
    files_sorter(path_to_256, append_prefix=True, design_labels=[],
                 faces=faces, fiducial_labels=[256], timesteps=1)

if run_moments:
    # Save moment arrays of the 256 cubes.

    moments_path = os.path.join(path_to_256, "moments")
    if not os.path.exists(moments_path):
        os.mkdir(moments_path)

    for face in fiducials_256:
        fid_name = fiducials_256[face][256][0]
        mask_mom = Mask_and_Moments(fid_name, scale=0.001 * u.K)
        mask_mom.make_moments()
        mask_mom.make_moment_errors()
        save_name = os.path.splitext(os.path.basename(fid_name))[0]

        mask_mom.to_fits(os.path.join(moments_256_path, save_name))


if run_distances:

    # Set which stats to run.
    statistics = copy(statistics_list)
    # statistics.remove("Dendrogram_Hist")
    # statistics.remove("Dendrogram_Num")
    statistics.append("DeltaVariance_Centroid_Curve")
    statistics.append("DeltaVariance_Centroid_Slope")

    statistics.remove("Tsallis")

    all_distances = {0: None, 1: None, 2: None}

    for face in fiducials.keys():
        print("On face {0} at {1}".format(face, time.ctime()))

        distances_storage = np.zeros((len(statistics), 5))

        if len(fiducials_256[face][256]) == 0:
            print("No face {} 256 fiducial found.".format(face))
            continue

        fid_256 = fiducials_256[face][256][0]
        dataset1 = \
            Mask_and_Moments.from_fits(fid_256,
                                       moments_path=moments_256_path).to_dict()

        fid_noise = 0.1 * np.nanpercentile(dataset1["cube"][0], 98)
        dendro_params_fid = {"min_value": 2 * fid_noise, "min_npix": 10}

        # Loop through 128 fiducials
        for i, fid_num in enumerate(ProgressBar(fiducials[face].keys())):
            fid_128 = fiducials[face][fid_num][-1]
            dataset2 = \
                Mask_and_Moments.from_fits(fid_128,
                                           moments_path=moments_path).to_dict()

            test_noise = 0.1 * np.nanpercentile(dataset2["cube"][0], 98)
            noise_value = [fid_noise, test_noise]

            dendro_params_test = {"min_value": 2 * test_noise, "min_npix": 80}
            dendro_params = [dendro_params_fid, dendro_params_test]

            # Pass the inertial ranges and spatial ranges for fitting
            # The 256 inertial range is k=5 to k~22
            if dataset1['moment0'][0].shape[0] == 256:
                low_range = [4.5 / 256., 4.5 / 128.]
                high_range = [22. / 256., 16. / 128.]
                low_spat_range = [11.] * 2
                high_spat_range = [51., 27.]
            else:
                low_range = [4.5 / 128.] * 2
                high_range = [16. / 128.] * 2
                low_spat_range = [11.] * 2
                high_spat_range = [27.] * 2

            inertial_range = [low_range, high_range]
            spatial_range = [low_spat_range, high_spat_range]

            print("Inertial Range: {}".format(inertial_range))
            print("Spatial Range: {}".format(spatial_range))

            vcs_break = -0.8
            vca_break = None

            # Downgrade the cubes to have 100 spectral channels. Small scales
            # are dominated by thermal motion, and so act like noise
            vcs_regrid = [100, 100]

            if i == 0:
                distances, fiducial_models = \
                    stats_wrapper(dataset1, dataset2,
                                  statistics=statistics,
                                  dendro_params=dendro_params,
                                  inertial_range=inertial_range,
                                  spatial_range=spatial_range,
                                  vcs_break=vcs_break,
                                  vca_break=vca_break,
                                  periodic_bounds=[True, True],
                                  vcs_regrid=vcs_regrid)
                all_fiducial_models = fiducial_models
            else:
                distances = \
                    stats_wrapper(dataset1, dataset2,
                                  fiducial_models=all_fiducial_models,
                                  multicore=True,  # don't need fiducial models
                                  statistics=statistics,
                                  dendro_params=dendro_params,
                                  inertial_range=inertial_range,
                                  spatial_range=spatial_range,
                                  vcs_break=vcs_break,
                                  vca_break=vca_break,
                                  periodic_bounds=[True, True],
                                  vcs_regrid=vcs_regrid)
            distances = [distances]
            distances_storage[:, i:i + 1] = \
                sort_distances(statistics, distances).T

        all_distances[face] = distances_storage

        # Save the results
        df = DataFrame(all_distances[face], index=statistics).T

        if run_on_regridded:
            output_name = \
                os.path.join(output_path,
                             "rescompare_regridded_face{}.csv".format(face))
        else:
            output_name = os.path.join(output_path,
                                       "rescompare_face{}.csv".format(face))

        df.to_csv(output_name)

if run_analysis:

    # This grabs the most recent clean analysis
    clean_path = os.path.join(results_path, "clean")
    folds = [os.path.join(clean_path, fold) for fold in
             os.listdir(clean_path) if "clean" in fold]
    folds.sort()
    print("Using the clean results from: {}".format(folds[-1]))
    path_to_128dists = folds[-1]


    niters = 10000

    for face in faces:
        print("Running on face {}.".format(face))

        if run_on_regridded:
            filename = \
                os.path.join(output_path,
                             "rescompare_regridded_face{}.csv".format(face))
        else:
            filename = os.path.join(output_path,
                                    "rescompare_face{}.csv".format(face))

        fids256 = \
            read_csv(filename, index_col=0)
        # Now load in the fiducial-fiducial distances of the 128 cubes
        fids128 = \
            read_csv(os.path.join(path_to_128dists,
                                  "fiducials_{0}_{0}.csv".format(face)),
                     index_col=0)[fids256.columns]

        pvals = dict.fromkeys(fids256.columns)

        for stat in ProgressBar(fids256.columns):
            # Stack all distances together w/ a dummy variable
            dists = np.hstack([fids256[stat], fids128[stat]])
            x = np.hstack([np.ones_like(fids256[stat]),
                           np.zeros_like(fids128[stat])])
            x = sm.add_constant(x)

            model = sm.OLS(dists, x)
            results = model.fit()
            mean_diff = results.params[1]

            # Now permute and re-test
            counter = 0
            for _ in xrange(niters):
                perm_dists = np.random.permutation(dists)
                perm_results = sm.OLS(perm_dists, x).fit()
                if perm_results.params[1] > mean_diff:
                    counter += 1
            pvals[stat] = float(counter) / float(niters)

        pval_df = DataFrame(pvals, index=[0])

        if run_on_regridded:
            output_name = \
                os.path.join(output_path,
                             "pvals_regridded_face_{}.csv".format(face))
        else:
            output_name = os.path.join(output_path,
                                       "pvals_face_{}.csv".format(face))
        pval_df.to_csv(output_name)
