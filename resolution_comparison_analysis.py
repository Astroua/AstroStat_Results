
from turbustat.data_reduction import Mask_and_Moments
from turbustat.statistics import stats_wrapper, statistics_list
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

from jasper.analysis_funcs import files_sorter, sort_distances

import matplotlib.pyplot as p
import seaborn as sns
sns.set_style("ticks")

p.ioff()

np.random.seed(446784788)

# Note that the moments need to be saved. See reduce_and_save_moments.py
# in jasper/

run_regrid = False
# To run on the regridded cubes, enable this
try:
    run_on_regridded = True if sys.argv[1] == "T" else False
except IndexError:
    run_on_regridded = True
run_distances = True
run_analysis = True

# path_to_data = "/media/eric/Data_3/Astrostat/SimSuite8/"
path_to_data = "/media/eric/Data_3/Astrostat/Fiducial_reproc/"
moments_path = os.path.join(path_to_data, "moments/")

# Only running on face 0
faces = [0]

if run_on_regridded:
    print("Running on regridded 256 cubes.")
    path_to_256 = "/media/eric/Data_3/Astrostat/Fiducial_256/regrid/"
else:
    print("Running on normal 256 cubes.")
    path_to_256 = "/media/eric/Data_3/Astrostat/Fiducial_256/"
moments_256_path = os.path.join(path_to_256, "moments/")

path_to_128dists = \
    "/media/eric/Data_3/Astrostat/results/clean/clean_20160725101641046250"

output_path = "/media/eric/Data_3/Astrostat/results/resolution_comparison"

if run_regrid:
    # Regrid the 256 cubes to 128 and save them.
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

        regrid_cube.write(os.path.join(path_to_256, "regrid",
                                       fid_256.split("/")[-1]))

if run_distances:

    # Sort those from the 128 set, and keep only the first timestep
    fiducials, _, _ = \
        files_sorter(path_to_data, timesteps=1, faces=faces,
                     append_prefix=True, design_labels=[])

    # Now the 256 cubes
    fiducials_256, _, _ = \
        files_sorter(path_to_256, append_prefix=True, design_labels=[],
                     faces=faces, fiducial_labels=[256], timesteps=1)

    # Set which stats to run.
    statistics = copy(statistics_list)
    # statistics.remove("Dendrogram_Hist")
    # statistics.remove("Dendrogram_Num")

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

        # Loop through 128 fiducials
        for i, fid_num in enumerate(fiducials[face].keys()):
            fid_128 = fiducials[face][fid_num][0]
            dataset2 = \
                Mask_and_Moments.from_fits(fid_128,
                                           moments_path=moments_path).to_dict()

            if i == 0:
                distances, fiducial_models = \
                    stats_wrapper(dataset1, dataset2,
                                  statistics=statistics)
                all_fiducial_models = fiducial_models
            else:
                distances = \
                    stats_wrapper(dataset1, dataset2,
                                  fiducial_models=all_fiducial_models,
                                  statistics=statistics)
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
