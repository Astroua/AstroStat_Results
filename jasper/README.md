# Scripts Description

This directory has the scripts used to run the distance comparisons for the set of simulated data cubes. The Jasper cluster used Torque as its job scheduler and the job submission scripts (`.pbs`) all use this format.

The scripts have the option to use MPI for running on multiple nodes. We don't recommend using this setting as it lead to various job crashes. This likely was due to issues with the Jasper cluster itself (it was near its defunding date when this analysis was run), but we can't guarantee the results using this option.

Below are descriptions of the vital scripts used.

* `reduce_and_save_moments.py` -- A function for masking the data cubes (primarily where noise has been added) and saving the masked cubes with their respective moment maps. Calls the normalizing/regridding function from `preprocessing.py`
* `make_moments.pbs`, `make_moments_noise.pbs`, `make_moments_regrid.pbs`, `mask_moments_complete.pbs` -- These submit jobs to pre-compute and save the masked cubes and their respective moment (and moment uncertainty) maps.
* `run_all.sh` -- This runs everything! The resource requirements and which analysis jobs to submit are controlled at the top of this script. Many of these analysis were run in two modes: an average over all time steps, and a comparison only at the common free-fall times. Both versions of the analysis will be submitted automatically. Note that the file structure is specifically set and is not generalized.
* `fiducial_submit.pbs` -- `run_all.sh` passes arguments to this job file, and submits it. It is essentially a job template that can be filled in.
* `output_mpi_load.py` -- This is the script called through run_all.sh. It's role is to organize the cube names to determine their time step and order in the experimental design. While the name includes "mpi", it can be disabled for running on a single cluster node. `output.py` and `output_mpi.py` are older versions with slight variations and should not be used.
* `preprocessor.py` -- Contains the function for normalizing and reprojecting the data cubes to have a common maximum intensity and average line width. See the Appendix section of the paper for a full description. Note that, since the spectral axis has been redefined from its original, the VCS is expected to have numerical effects from the resulting data cubes.
* `complete_comparisons.py` -- Like `output_mpi_load.py`, but for comparisons with the COMPLETE data.
* `complete_to_complete.pbs`, `complete_to_fid.pbs`, `des_to_complete.pbs` -- Job files for the COMPLETE comparisons. Each uses `complete_comparisons.py`.
* `run_viewing_angle_distance.py`, `run_viewing_angle_distance.pbs` -- Comparisons between a simulated observed cube from two different directions. This was to determine the effect of the viewing angle on the distance metrics. **NOT USED IN THIS PAPER**.
* `wrapping_function.py` -- Contains a function that runs all of the TurbuStat statistics on two different data cubes (and their moments). Their are a number of added inputs to account for variations in the analysis. Key to this is the handling of boundaries in the data (periodic?), the noise levels, loading of pre-saved results, and the inertial ranges to be used when fitting power-law indices.