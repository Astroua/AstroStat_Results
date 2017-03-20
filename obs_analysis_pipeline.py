
'''
 Reproduce the sim vs obs plots from the paper

Expects the following folder structure before running:
    - Des_to_Obs
        - HDF5
    - Obs_to_Fid
        - HDF5
    - Obs_to_Obs

The design matrix (saved in csv format) should be in the parent directory
(given as sys.argv[1]).

The obs to fiducial comparisons require the noise added simulation results
to be copied. The path to these results is given by sys.argv[2].

'''

# Create data tables of the observational results
# Run from Dropbox/AstroStatistics/Full Factorial/Observational Results/

import os
import shutil
import sys

from analysis.convert_results import concat_convert_HDF5
from analysis import convert_format, comparison_plot


path_to_data = os.path.abspath(sys.argv[1])
path_to_noise_sim_results = os.path.abspath(sys.argv[2])
design_matrix = os.path.join("/".join(path_to_data.split("/")[:-1]),
                             "Design7Matrix.csv")
# Pass as 0,2 for multiple faces
faces = [int(face) for face in sys.argv[3].split(",")]

# Obs to Fids
obsfid_path = os.path.join(path_to_data, "Obs_to_Fid")

# Des to Obs
desobs_path = os.path.join(path_to_data, "Des_to_Obs")

for face in faces:

    # Obs to Fids
    convert_format(os.path.join(obsfid_path, "HDF5"), face)
    shutil.move(os.path.join(obsfid_path,
                             "HDF5/complete_distances_face_{}.csv".format(face)),
                obsfid_path)

    # Now move the noisy simulation comparisons into the folder
    shutil.copy(os.path.join(path_to_noise_sim_results,
                             "distances_{0}_{0}.csv".format(face)),
                obsfid_path)
    shutil.copy(os.path.join(path_to_noise_sim_results,
                             "fiducials_{0}_{0}.csv".format(face)),
                obsfid_path)

    # Des to Obs
    concat_convert_HDF5(os.path.join(desobs_path, "HDF5"), face=face,
                        interweave=True, average_axis=0)
    shutil.move(os.path.join(desobs_path,
                             "HDF5/distances_{}.csv".format(face)),
                desobs_path)
    shutil.move(os.path.join(desobs_path, "distances_{}.csv".format(face)),
                os.path.join(desobs_path, "distances_{}_obs.csv".format(face)))

    # Finished pipeline should have a similar conversion for Obs_to_Obs, but
    # this remains manual for now. These are copied into Des_to_Obs, where the
    # are the pseudo-fiducials

    # Same "fiducials" for both face comparisons
    obsobs_path = os.path.join(path_to_data, "Obs_to_Obs")
    shutil.copy(os.path.join(obsobs_path, "complete_comparisons.csv"),
                os.path.join(desobs_path, "fiducial_{}_obs.csv".format(face)))

# Now make the plots

comparisons = ["{0}_{0}".format(face) for face in faces]
comparisons_obs = ["{0}_obs".format(face) for face in faces]

# Obs to Fid

if not os.path.exists(os.path.join(obsfid_path, "Distance_Plots")):
    os.mkdir(os.path.join(obsfid_path, "Distance_Plots"))

comparison_plot(
    obsfid_path, comparisons=comparisons,
    out_path=os.path.join(obsfid_path, "Distance_Plots"),
    num_fids=5, design_matrix=design_matrix, obs_to_fid=True, legend=False,
    obs_legend=True,
    statistics=["Cramer", "DeltaVariance_Curve", "DeltaVariance_Slope",
                "DeltaVariance_Centroid_Curve", "DeltaVariance_Centroid_Slope",
                "Dendrogram_Hist", "Dendrogram_Num", "PCA", "SCF", "VCA",
                "VCS", "VCS_Small_Scale", "VCS_Large_Scale",
                "Skewness", "Kurtosis"],
    show_title=True if len(faces) == 1 else False,
    use_tightlayout=True)

# Des to Obs

if not os.path.exists(os.path.join(desobs_path, "Distance_Plots")):
    os.mkdir(os.path.join(desobs_path, "Distance_Plots"))

comparison_plot(
    desobs_path, comparisons=comparisons_obs,
    out_path=os.path.join(desobs_path, "Distance_Plots"),
    num_fids=3, design_matrix=design_matrix,
    legend_labels=["Oph A", "IC 348", "NGC 1333"],
    statistics=["Cramer", "DeltaVariance_Curve", "DeltaVariance_Slope",
                "DeltaVariance_Centroid_Curve", "DeltaVariance_Centroid_Slope",
                "Dendrogram_Hist", "Dendrogram_Num", "PCA", "SCF", "VCA",
                "VCS", "VCS_Small_Scale", "VCS_Large_Scale",
                "Skewness", "Kurtosis"],
    show_title=True if len(faces) == 1 else False,
    use_tightlayout=True)

