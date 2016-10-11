
import sys
from glob import glob
import os
from pandas import DataFrame, read_csv
import statsmodels.api as sm
from astropy.utils.console import ProgressBar
import numpy as np


try:
    path_to_data = sys.argv[1]
except IndexError:
    path_to_data = raw_input("Provide path to face comparison data: ")

save = True
niters = 5000

data_files = np.array(glob(os.path.join(path_to_data, "view_angle_*.csv")))

# There should be 3 comparisons
assert len(data_files) == 3

which_0_1 = np.array([True if "0_1" in file else False for file in data_files])
face_0_1 = read_csv(data_files[which_0_1])

which_0_2 = np.array([True if "0_2" in file else False for file in data_files])
face_0_2 = read_csv(data_files[which_0_2])

which_1_2 = np.array([True if "1_2" in file else False for file in data_files])
face_1_2 = read_csv(data_files[which_1_2])

comparisons = [[face_0_1, face_0_2, "0_1-0_2"],
               [face_0_1, face_1_2, "0_1-1_2"],
               [face_0_2, face_1_2, "0_2-1_2"]]

face_comp_labels = ["0_1-0_2", "0_1-1_2", "0_2-1_2"]
pvals = dict.fromkeys(face_comp_labels)

stats = face_0_1.columns.copy()
stats = stats.drop(["Data 1", "Data 2"])


for data1, data2, label in comparisons:

    pvals[label] = {}

    for stat in ProgressBar(stats):
        # Stack all distances together w/ a dummy variable
        dists = np.hstack([data1[stat], data2[stat]])
        x = np.hstack([np.ones_like(data1[stat]),
                       np.zeros_like(data2[stat])])
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

pval_df = DataFrame(pvals, index=face_comp_labels)

if save:
    pval_df.T.to_csv(os.path.join(path_to_data, "view_angle_pvals.csv"))

# Look for which comparisons show significant distance properties with
# different faces
crit_p = 0.05
signif_df = pval_df < crit_p

if save:
    signif_df.T.to_csv(os.path.join(path_to_data,
                                    "view_angle_signifdiff.csv"))

# try:
#     path_to_data = sys.argv[1]
# except IndexError:
#     path_to_data = raw_input("Provide path to fiducial data: ")

# fiducials = glob(os.path.join(path_to_data, "fiducials*.csv"))

# if not len(fiducials) == 9:
#     raise Warning("All face comparisons require there be 9 fiducial files. "
#                   "Only {} were found".format(len(fiducials)))

# faces = ["0_0", "0_1", "0_2", "1_0", "1_1", "1_2", "2_0", "2_1", "2_2"]
# fiducial_data = dict.fromkeys(faces)

# for file in fiducials:
#     if not os.path.basename(file)[-7:-4] in faces:
#         raise Warning("{} does not match face combination pattern."
#                       .format(file))

#     for face in faces:
#         if face in os.path.basename(file):
#             fiducial_data[face] = read_csv(file, index_col=0)


# # Comparisons between to run p-tests on.
# comparisons = [["0_1", "1_0"], ["0_2", "2_0"], ["1_2", "2_1"], ["0_0", "1_1"],
#                ["0_0", "2_2"], ["2_2", "1_1"]]

# all_pvals = []

# save = True
# niters = 5000

# # We're going to assume that each has the same columns. This should be true
# # for a complete set of comparisons.
# stats = fiducial_data["0_0"].columns.copy()
# stats = stats.drop(["Fiducial 1", "Fiducial 2"])
# pvals = dict.fromkeys(stats)
# for key in pvals:
#     pvals[key] = np.empty((6,))

# # Track the comparisons
# pvals["Face Pair 1"] = np.empty((6,), dtype=object)
# pvals["Face Pair 2"] = np.empty((6,), dtype=object)

# for i, comp in enumerate(comparisons):

#     data1 = fiducial_data[comp[0]]
#     data2 = fiducial_data[comp[1]]

#     for stat in ProgressBar(stats):
#         # Stack all distances together w/ a dummy variable
#         dists = np.hstack([data1[stat], data2[stat]])
#         x = np.hstack([np.ones_like(data1[stat]),
#                        np.zeros_like(data2[stat])])
#         x = sm.add_constant(x)

#         model = sm.OLS(dists, x)
#         results = model.fit()
#         mean_diff = results.params[1]

#         # Now permute and re-test
#         counter = 0
#         for _ in xrange(niters):
#             perm_dists = np.random.permutation(dists)
#             perm_results = sm.OLS(perm_dists, x).fit()
#             if perm_results.params[1] > mean_diff:
#                 counter += 1
#         pvals[stat][i] = float(counter) / float(niters)

#     pvals["Face Pair 1"][i] = comp[0]
#     pvals["Face Pair 2"][i] = comp[1]

# pval_df = DataFrame(pvals, index=np.arange(6))

# if save:
#     pval_df.T.to_csv(os.path.join(path_to_data, "fiducial_view_pvals.csv"))

# # Now let's try to reduce this to "viewing angle alone" vs.
# # "magnetic field direction"
# crit_p = 0.1
# sensitivity = dict.fromkeys(stats)
# for key in sensitivity:
#     sensitivity[key] = np.empty((2,))

#     # Viewing angle
#     # Is there a difference between 0_1 and 1_0 [0]
#     # Also 0_0 or 1_1? [3]
#     # sensitivity[key][0] = np.logical_and(pvals[key][0] < crit_p,
#     #                                      pvals[key][3] < crit_p)
#     sensitivity[key][0] = pvals[key][0] < crit_p

#     # Possible difference from B-field direction
#     # Is there a difference b/w 0_2, 2_0 [1] and 1_2, 2_1 [2].
#     # And maybe 0_0, 2_2 [4] and 1_1, 2_2 [5]?
#     first_pair = np.logical_and(pvals[key][1] < crit_p,
#                                 pvals[key][2] < crit_p)
#     second_pair = np.logical_and(pvals[key][4] < crit_p,
#                                  pvals[key][5] < crit_p)

#     # sensitivity[key][1] = np.logical_and(first_pair, second_pair)
#     sensitivity[key][1] = first_pair

# sensitivity_df = DataFrame(sensitivity,
#                            index=["Viewing Angle", "B-field Direction"])

# if save:
#     sensitivity_df.T.to_csv(os.path.join(path_to_data,
#                                          "fiducial_view_sensitivity.csv"))
