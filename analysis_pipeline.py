
import numpy as np
import analysis as ta
from turbustat.statistics import statistics_list
import os
import subprocess
import sys
import shutil
from datetime import datetime

import seaborn as sb

# Set ticks facing inward
sb.set_context("paper")
sb.set(font='Times New Roman', style='ticks', font_scale=1.0)
sb.set_style({"xtick.direction": "in","ytick.direction": "in"})

'''
Runs the basic analysis pipeline, starting from the outputted HDF5 files

Creates a time-stamped folder with the expected structure.
'''


def timestring():
    TimeString = datetime.now().strftime("%Y%m%d%H%M%S%f")
    return TimeString


# Specify path with results
path = os.path.abspath(sys.argv[1])

if path.split("/")[-1] == "HDF5_files":
    hdf5_path = path + "/"
    path = "/".join(path.split("/")[:-1]) + "/"
    if path == "/":
        path = "./"
else:
    hdf5_path = os.path.join(path, "HDF5_files/")
    path += "/"

design_matrix = sys.argv[2]

if design_matrix == "None":
    design_matrix = None

scripts_path = sys.argv[3]

run_name = sys.argv[4]

tstep_choice = sys.argv[5]

added_noise = True if sys.argv[6] == "T" else False

faces = sys.argv[7].split()

faces = ["{0}_{0}".format(face) for face in faces]

statistics = statistics_list
statistics.append("DeltaVariance_Centroid_Slope")
statistics.append("DeltaVariance_Centroid_Curve")

if tstep_choice == "mean":
    mode = 'mean'
    des_tsteps = None
    fid_tsteps = None
# elif tstep_choice == "freefall":
#     mode = 'choice'
#     # Timesteps for Simulation set 8 (used in the Tools 2016 paper).
#     des_tsteps = \
#         np.array([25, 26, 21, 21, 23, 25, 21, 21, 25, 26, 21, 22, 23, 26, 21,
#                   21, 30, 30, 27, 28, 30, 30, 27, 23, 30, 30, 27, 25, 30, 30,
#                   24, 25])
#     # Make sure this matches the number of design sims
#     assert len(des_tsteps) == 32

#     fid_tsteps = np.array([25, 24, 26, 26, 26])
#     assert len(fid_tsteps) == 5

#     # Subtract the 'zeroth' timestep to get the right index
#     des_tsteps -= 21
#     fid_tsteps -= 21

else:
    # raise ValueError("tstep_choice must be 'mean' or 'freefall'")
    raise ValueError("tstep_choice must be 'mean'.")

# Redefine the path with the timestamped output folder
path = os.path.join(path, run_name + "_" + timestring())
os.mkdir(path)
shutil.copytree(hdf5_path, os.path.join(path, "HDF5_files"))
hdf5_path = os.path.join(path, "HDF5_files")

os.chdir(path)

# Convert into combined csv files.

good_comparison = []

print("Converting to combined csv files.")

for face1 in [0, 1, 2]:
    for face2 in [0, 1, 2]:
        try:
            ta.convert_format(hdf5_path, face1, face2=face2,
                              design=design_matrix,
                              mode=mode, tsteps=des_tsteps)
            shutil.move(os.path.join(hdf5_path,
                                     "distances_{0}_{1}.csv".format(face1,
                                                                    face2)),
                        path)
            good_comparison.append("{0}_{1}".format(face1, face2))
        except StandardError as err:
            print err

# Next convert the fiducial comparisons

for fil in os.listdir(hdf5_path):
    full_file = os.path.join(hdf5_path, fil)
    if os.path.isfile(full_file) and "fid_comp" in fil:
        print("Fiducial comparison: " + full_file)
        out_name = ta.convert_fiducial(full_file,
                                       return_name=True,
                                       mode=mode,
                                       tsteps=fid_tsteps)

        # shutil.move(out_name, path)

# Now make the distance plots.

print("Making distance plots.")

# Check that the requested face files were found
for face in faces:
    if face in good_comparison:
        continue
    raise ValueError("Face {0} not found. Those found are:"
                     " {1}".format(face, good_comparison))

if not os.path.exists(os.path.join(path, "Distance_Plots")):
    os.mkdir(os.path.join(path, "Distance_Plots"))

ta.comparison_plot(path, comparisons=faces,
                   out_path=os.path.join(path, "Distance_Plots"),
                   design_matrix=design_matrix,
                   statistics=statistics,
                   show_title=True if len(faces) == 1 else False,
                   use_tightlayout=True)

if not os.path.exists(os.path.join(path, "Distance_Plots_Paper")):
    os.mkdir(os.path.join(path, "Distance_Plots_Paper"))

ta.comparison_plot(path, comparisons=faces,
                   out_path=os.path.join(path, "Distance_Plots_Paper"),
                   design_matrix=design_matrix,
                   num_fids=5,
                   statistics=statistics,
                   show_title=True if len(faces) == 1 else False,
                   use_tightlayout=True)

# Run the R-script to fit the data to the model

os.chdir(path)

# This isn't robust. Assumes that when only 1 face is used, it is the 0 0
# comparison
if len(faces) == 1:
    print("Only found one face: {}".format(faces))
    factorial_script = "FactorialAnalysis_oneface.R"
    noise_script = "noise_validation_oneface.r"
    signal_script = "signal_validation_oneface.r"
    params = ["pb", "m", "k", "sf", "vp"]
else:
    print("Found multiple faces: {}".format(faces))
    factorial_script = "FactorialAnalysis.R"
    noise_script = "noise_validation.r"
    signal_script = "signal_validation.r"
    params = ["fc", "pb", "m", "k", "sf", "vp"]


print("Fitting model of given design.")

subprocess.call(['Rscript',
                 os.path.join(scripts_path, factorial_script)])

# This should create two output tables of the whole dataset.

# Now run the metric validation

print("Running metric validation.")

subprocess.call(['Rscript',
                 os.path.join(scripts_path, noise_script),
                 "10000"])

subprocess.call(['Rscript',
                 os.path.join(scripts_path, signal_script),
                 "10000"])

# Finally, create the model plots

print("Creating model plots.")

if not os.path.exists(os.path.join(path, "Model_Plots")):
    os.mkdir(os.path.join(path, "Model_Plots"))

# What are we considering the min t-value to be?
# min_tvalue = 2.0  # 95%
min_tvalue = 3.46  # 99.9%

ta.effect_plots("DataforFits.csv", "ResultsFactorial.csv", save=True,
                out_path='Model_Plots', min_tvalue=min_tvalue,
                params=params)

import seaborn as sb
sb.set_context("paper")
sb.set(font='Times New Roman', style='ticks', font_scale=1.2)

figsize = (4.2, 6.0)

# Coef plots for all terms
ta.make_coefplots("DataforFits.csv", save=True, out_path="Model_Plots",
                  min_tvalue=min_tvalue, statistics=statistics,
                  endog_formula="*".join(params),
                  figsize=figsize)

# Only show results of the good statistics
if added_noise:
    good_stats = ["Cramer", "DeltaVariance_Centroid_Curve",
                  "DeltaVariance_Centroid_Slope", "PDF_Lognormal",
                  "SCF", "VCA", "VCS_Large_Scale",
                  "Skewness", "Kurtosis",
                  "Genus", "Bispectrum", "PSpec",
                  "MVC", "PCA", "Dendrogram_Num"]
else:
    good_stats = ["Cramer", "DeltaVariance_Centroid_Curve",
                  "DeltaVariance_Centroid_Slope", "PDF_Lognormal",
                  "SCF", "VCA", "VCS", "VCS_Small_Scale", "VCS_Large_Scale",
                  "Skewness", "Kurtosis",
                  "Genus", "Bispectrum", "PSpec",
                  "MVC", "PCA", "Dendrogram_Num"]

sb.set_context("paper")
sb.set(font='Times New Roman', style='ticks', font_scale=1.0)
sb.set_style({"xtick.direction": "out", "ytick.direction": "out"})

# Set the figure size to be 1 column w/ height around half the page
figsize = (4.2, 6.0)

ta.map_all_results("ResultsFactorial.csv", normed=False, max_order=None,
                   save_name="map_all_results.pdf",
                   out_path='Model_Plots', statistics=good_stats,
                   min_tvalue=min_tvalue, max_tvalue=30,
                   figsize=figsize)

print("Finished!")
