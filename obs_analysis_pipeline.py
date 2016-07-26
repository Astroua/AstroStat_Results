
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
design_matrix = os.path.join("/".join(path_to_data.split("/")[:-1]), "Design7Matrix.csv")

# Obs to Fids

obsfid_path = os.path.join(path_to_data, "Obs_to_Fid")

# Face 0
convert_format(os.path.join(obsfid_path, "HDF5"), 0)
shutil.move(os.path.join(obsfid_path, "HDF5/complete_distances_face_0.csv"),
            obsfid_path)

# Face 1
# convert_format(hdf5_path, 1)
# shutil.move(hdf5_path+"complete_distances_face_1.csv", path)

# Face 2
convert_format(os.path.join(obsfid_path, "HDF5"), 2)
shutil.move(os.path.join(obsfid_path, "HDF5/complete_distances_face_2.csv"),
            obsfid_path)

# Now move the noisy simulation comparisons into the folder
shutil.copy(os.path.join(path_to_noise_sim_results, "distances_0_0.csv"),
            obsfid_path)
shutil.copy(os.path.join(path_to_noise_sim_results, "distances_2_2.csv"),
            obsfid_path)
shutil.copy(os.path.join(path_to_noise_sim_results, "fiducials_0_0.csv"),
            obsfid_path)
shutil.copy(os.path.join(path_to_noise_sim_results, "fiducials_2_2.csv"),
            obsfid_path)

# Des to Obs

desobs_path = os.path.join(path_to_data, "Des_to_Obs")

# Face 0
concat_convert_HDF5(os.path.join(desobs_path, "HDF5"), face=0,
                    interweave=True, average_axis=0)
shutil.move(os.path.join(desobs_path, "HDF5/distances_0.csv"), desobs_path)
shutil.move(os.path.join(desobs_path, "distances_0.csv"),
            os.path.join(desobs_path, "distances_0_obs.csv"))

# Face 2
concat_convert_HDF5(os.path.join(desobs_path, "HDF5"), face=2,
                    interweave=True, average_axis=0)
shutil.move(os.path.join(desobs_path, "HDF5/distances_2.csv"), desobs_path)
shutil.move(os.path.join(desobs_path, "distances_2.csv"),
            os.path.join(desobs_path, "distances_2_obs.csv"))

# Finished pipeline should have a similar conversion for Obs_to_Obs, but this
# remains manual for now. These are copied into Des_to_Obs, where the are the
# pseudo-fiducials

# Same "fiducials" for both face comparisons
obsobs_path = os.path.join(path_to_data, "Obs_to_Obs")
shutil.copy(os.path.join(obsobs_path, "complete_comparisons.csv"),
            os.path.join(desobs_path, "fiducial_0_obs.csv"))
shutil.copy(os.path.join(obsobs_path, "complete_comparisons.csv"),
            os.path.join(desobs_path, "fiducial_2_obs.csv"))

# Now make the plots

# Obs to Fid

obs_to_fid_comparisons = ["0_0", "2_2"]

if not os.path.exists(os.path.join(obsfid_path, "Distance_Plots")):
    os.mkdir(os.path.join(obsfid_path, "Distance_Plots"))

comparison_plot(
    obsfid_path, comparisons=obs_to_fid_comparisons,
    out_path=os.path.join(obsfid_path, "Distance_Plots"),
    num_fids=5, design_matrix=design_matrix, obs_to_fid=True, legend=False,
    obs_legend=True,
    statistics=["Cramer", "DeltaVariance", "Dendrogram_Hist",
                "Dendrogram_Num", "PCA", "SCF", "VCA", "VCS",
                "VCS_Small_Scale", "VCS_Large_Scale", "VCS_Break",
                "Skewness", "Kurtosis"])

# Des to Obs

if not os.path.exists(os.path.join(desobs_path, "Distance_Plots")):
    os.mkdir(os.path.join(desobs_path, "Distance_Plots"))

des_to_obs_comparisons = ["0_obs", "2_obs"]

comparison_plot(
    desobs_path, comparisons=des_to_obs_comparisons,
    out_path=os.path.join(desobs_path, "Distance_Plots"),
    num_fids=3, design_matrix=design_matrix,
    legend_labels=["Oph A", "IC 348", "NGC 1333"],
    statistics=["Cramer", "DeltaVariance", "Dendrogram_Hist",
                "Dendrogram_Num", "PCA", "SCF", "VCA", "VCS",
                "VCS_Small_Scale", "VCS_Large_Scale", "VCS_Break",
                "Skewness", "Kurtosis"])
