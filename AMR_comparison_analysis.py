
'''
Compare the timestep 30 fiducials with and without AMR.
'''

import numpy as np
import astropy.units as u
import os
import sys
from astropy.utils.console import ProgressBar
import pandas as pd
from copy import copy
import statsmodels.api as sm
import statsmodels.formula.api as smf

from turbustat.data_reduction import Mask_and_Moments
from turbustat.statistics import statistics_list

from jasper.analysis_funcs import (files_sorter, sort_distances,
                                   timestep_wrapper)


# sys.argv[1]  - path to the full (non-AMR) dataset
# sys.argv[2] - path to the AMR fiducials
# sys.argv[3] - output path to save results in

path_to_data = sys.argv[1]
moments_path = os.path.join(path_to_data, "moments/")

path_to_amrdata = sys.argv[2]
amrmoments_path = os.path.join(path_to_amrdata, "moments/")

results_path = sys.argv[3]
save_path = os.path.join(results_path, "amr_comparison")
if not os.path.exists(save_path):
    os.mkdir(save_path)

faces = [0] # [0, 2]

# fiducials, _, _ = \
#     files_sorter(path_to_data, timesteps='last', faces=faces,
#                  append_prefix=True, design_labels=[], verbose=False)

# # Now the AMR cubes
# fiducials_amr, _, _ = \
#     files_sorter(path_to_amrdata, append_prefix=True, design_labels=[],
#                  faces=faces, timesteps='last', verbose=False)

# # If the AMR moments path doesn't exist, make the moment arrays and save.
# if not os.path.exists(amrmoments_path):

#     os.mkdir(amrmoments_path)

#     for face in faces:
#         for fid in fiducials_amr[face]:
#             fid_name = fiducials_amr[face][fid]
#             mask_mom = Mask_and_Moments(fid_name,
#                                         scale=0.001 * u.K)
#             mask_mom.make_moments()
#             mask_mom.make_moment_errors()

#             save_name = os.path.splitext(os.path.basename(fid_name))[0]
#             mask_mom.to_fits(os.path.join(amrmoments_path, save_name))


# Now run the distances AMR vs. none.
statistics = copy(statistics_list)
statistics.append("DeltaVariance_Centroid_Curve")
statistics.append("DeltaVariance_Centroid_Slope")

print "Statistics to run: %s" % (statistics)
num_statistics = len(statistics)

# for face in faces:
#     for fid in ProgressBar([3, 4]):

#         distances = np.zeros((len(statistics),
#                               len(fiducials_amr[face].keys())))

#         fid_name = fiducials[face][fid]

#         for i, amr_fid in enumerate(fiducials_amr[face].keys()):
#             fid_amr_name = fiducials_amr[face][amr_fid]

#             out = timestep_wrapper(fid_name, fid_amr_name, statistics, False)

#             out = [out]
#             distances[:, i] = sort_distances(statistics, out).T.squeeze()

#         df = pd.DataFrame(distances, index=statistics).T

#         # Save each fiducial as a csv file
#         save_name = "SimSuite8_fiducial{0}_amr_comparison_" \
#             "face_{1}.csv".format(fid, face)
#         df.to_csv(os.path.join(save_path, save_name))

#     # And the distances between the AMR fiducials.

#     dists = {}

#     for fid in ProgressBar(fiducials_amr[face].keys()):

#         fid_name = fiducials_amr[face][fid]

#         for i, amr_fid in enumerate(fiducials_amr[face].keys()):
#             if amr_fid == fid:
#                 continue

#             fid_amr_name = fiducials_amr[face][amr_fid]

#             out = timestep_wrapper(fid_name, fid_amr_name, statistics, False)

#             out = [out]
#             dists["{0}_{1}".format(fid, amr_fid)] = \
#                 sort_distances(statistics, out).T.squeeze()

#     df_fids = pd.DataFrame(dists, index=statistics).T

#     # Save each fiducial as a csv file
#     save_name = "SimSuite8_fiducial_to_fiducial_amr_comparison_" \
#         "face_{0}.csv".format(face)
#     df_fids.to_csv(os.path.join(save_path, save_name))

# Now run the analysis.

# Need the path to the most recent clean run
clean_path = os.path.join(results_path, "clean")
folds = [os.path.join(clean_path, fold) for fold in
         os.listdir(clean_path) if "clean" in fold]
folds.sort()
print("Using the clean results from: {}".format(folds[-1]))
path_to_128dists = folds[-1]

niters = 10000

for face in faces:
    print("Running on face {}.".format(face))

    # non-AMR to non-AMR
    fid_fid = pd.read_csv(os.path.join(folds[-1],
                                       "fiducials_{0}_{0}.csv".format(face)))

    # AMR to AMR
    filename = "SimSuite8_fiducial_to_fiducial_amr_comparison_" \
        "face_{}.csv".format(face)
    amr_amr = pd.read_csv(os.path.join(save_path, filename))

    # Load in each fiducial comparison and stack
    for fid in np.arange(0, 5):

        name = "SimSuite8_fiducial{0}_amr_comparison_face_" \
            "{1}.csv".format(fid, face)
        df = pd.read_csv(os.path.join(save_path, name))
        # Add in row to identify which fiducial this is
        df["Fiducial"] = pd.Series(np.array([fid] * 5))

        if fid == 0:
            amr_fid = df
        else:
            amr_fid = pd.concat([amr_fid, df])

    # Signficance of difference between AMR to Fid and Fid to Fid
    pvals_amr_to_fid = dict.fromkeys(amr_amr.columns)
    # Signficance of difference between AMR to AMR and Fid to Fid
    # So does the AMR just plain have more scatter associated?
    pvals_fid_to_fid = dict.fromkeys(amr_amr.columns)

    for stat in ProgressBar(pvals_amr_to_fid.keys()):

        # Stack all distances together w/ a dummy variable
        # QUESTION: Treat the testing against the same as one group?
        dists = np.hstack([amr_fid[stat], amr_amr[stat],
                           fid_fid[stat]]).astype(float)
        x = np.hstack([2 * np.ones_like(amr_fid[stat], dtype=np.float),
                       np.ones_like(amr_amr[stat], dtype=np.float),
                       np.zeros_like(fid_fid[stat], dtype=np.float)])

        stacked = np.hstack([dists[:, np.newaxis], x[:, np.newaxis]])
        df_dists = pd.DataFrame(stacked, columns=['Dists', 'Categ'])

        model = smf.ols("Dists ~ C(Categ)", data=df_dists)
        results = model.fit()
        mean_diff_amrfid = results.params[1]
        mean_diff_fidfid = results.params[2]

        # Now permute and re-test
        counter_amrfid = 0
        counter_fidfid = 0
        for _ in xrange(niters):
            perm_dists = np.random.permutation(dists)
            stacked_perm = np.hstack([perm_dists[:, np.newaxis],
                                      x[:, np.newaxis]])
            df_perm_dists = \
                pd.DataFrame(stacked_perm, columns=['Dists', 'Categ'])

            perm_model = smf.ols("Dists ~ C(Categ)", data=df_perm_dists)
            perm_results = perm_model.fit()
            if perm_results.params[1] > mean_diff_amrfid:
                counter_amrfid += 1
            if perm_results.params[2] > mean_diff_fidfid:
                counter_fidfid += 1

        pvals_amr_to_fid[stat] = float(counter_amrfid) / float(niters)
        pvals_fid_to_fid[stat] = float(counter_fidfid) / float(niters)

    pval_df = pd.DataFrame({"pvals AMR to Fid": pvals_amr_to_fid,
                            "pvals Fid to Fid": pvals_fid_to_fid})

    output_name = os.path.join(save_path,
                               "pvals_face_{}.csv".format(face))
    pval_df.to_csv(output_name)
