
import numpy as np
import analysis as ta
from turbustat.statistics import statistics_list
import os
import subprocess
import sys
import shutil
from datetime import datetime

'''
Runs the analysis for comparisons between the hot fiducials,
starting from the outputted HDF5 files

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

scripts_path = sys.argv[2]

run_name = sys.argv[3]

tstep_choice = sys.argv[4]

added_noise = True if sys.argv[5] == "T" else False

faces = sys.argv[6].split()

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
                              design=None,
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
                   design_matrix=None,
                   statistics=statistics)

if not os.path.exists(os.path.join(path, "Distance_Plots_Paper")):
    os.mkdir(os.path.join(path, "Distance_Plots_Paper"))


ta.comparison_plot(path, comparisons=faces,
                   out_path=os.path.join(path, "Distance_Plots_Paper"),
                   design_matrix=None,
                   num_fids=5,
                   statistics=statistics,
                   show_title=True if len(faces) == 1 else False,
                   use_tightlayout=True)


os.chdir(path)

# This isn't robust. Assumes that when only 1 face is used, it is the 0 0
# comparison
if len(faces) == 1:
    print("Only found one face: {}".format(faces))
    noise_script = "noise_validation_oneface.r"
else:
    print("Found multiple faces: {}".format(faces))
    noise_script = "noise_validation.r"

# Now run the metric validation

print "Running metric validation."

subprocess.call(['Rscript',
                 os.path.join(scripts_path, noise_script),
                 "10000"])

print("Finished!")
