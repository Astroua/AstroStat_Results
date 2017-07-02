# AstroStat_Results
Scripts for reproducing the results of Koch et al. (2017) (ADD LINK)

The complete set of data cubes is available at: (ADD LINK). The results used in the analysis are available in csv format at: (ADD LINK).

This repository contains the scripts for reproducing all aspects of the analysis. The scripts used for the distance metric comparisons are contained in the `Jasper` folder. The analysis scripts for the distance comparisons are described below.

The scripts are best used as part of the pipeline. Most scripts will only work when the expected file-naming convention is maintained. 

# Requirements:

The analysis was run using python 2.7. The required python packages are:

* numpy
* scipy
* matplotlib
* astropy>1.0
* pandas
* statsmodels
* patsy
* astrodendro

The modelling of the distance metric results requires the R language. Your path needs to include `Rscript`, which is used to run the fitting from within the python scripts. A few R packages are also needed:

* xtable
* MASS
* lme4

# Scripts Description:
* Jasper scripts -- Main comparisons are run on the Jasper cluster. This includes time-averaged and freefall comparisons of the clean cubes, noise-added cubes, hot fiducials, and observational comparisons. All jobs can be submitted through `jasper/run_all.sh`.
* Smaller analyses -- We do a few more comparisons that are smaller in scope, and so are not setup for running on the cluster. This includes: the effect of AMR (`AMR_comparison_analysis.py`), and the resolution comparison (`resolution_comparison_analysis.py`; we do 2 runs with the 256 grid cubes, and another regridding back down to 128 pixels).
* Analysis -- The main analysis is run through `run_all_analysis.sh`. This script takes the output from the jobs run on the Jasper cluster, creates a file directory structure, then runs the analyses. The analyses to run are depend on the variables set at the beginning of the script and the paths to the data. The main analysis is run through `analysis_pipeline.py`. This will convert the HDF5 files into the csv files expected for use with the R fitting scripts.

    * AMR_comparison_analysis.py -- Calculates p-values of the fiducial comparisons with and without AMR post-processing enabled in RADMC3D
    * FactorialAnalysis.R -- Fits the experimental design matrix to the distance values to find the parameter sensitivity. Includes the direction through the simulation that RADMC3D was run in (the "face") as an added parameter.
    * FactorialAnalysis_oneface.R -- Same as above, without treatment of the face as a variable.
    * analysis_pipeline.py -- The analysis pipeline used for the time-averaged and free-fall time only comparisons. If trying to use these scripts on another dataset, start by examining this script.
    * hot_compare_analysis.py -- Calculates p-values of the fiducials at 10 K and those at 40 K. This is a measure of the effect of temperature on the behaviour of the statistics.
    * model_basic_properties.py -- Run a sensitivity analysis (similar to that in FactorialAnalysis.R) between the peak intensity, total intensity, and mean line width. This determines the correlations inherent in the experimental design.
    * noise_validation.r -- Calculates p-values by permuting Fiducial-to-Design and Fiducial-to-Fiducial distances. This calculates whether there is any signal above the noise in the Fiducial-to-Design distances.
    * noise_validation_oneface.r -- Same as above, without treatment of the face as a variable.
    * obs_analysis_pipeline.py -- Creates the comparison plots between the distances from the simulated observation set and the distances calculated with respect to the COMPLETE 13CO data (for NGC 1333, Oph A, and IC 348).
    * resolution_comparison_analysis.py -- Calculates p-values for the distances between the 256-pixel fiducials and the 128-pixel fiducials. This calculates whether resolution is a significant effect for any of the statistics.
    * run_all_analysis.sh -- The entire analysis can be run from this script. 
    * signal_validation.r -- Calculates p-values by permuting the Fiducial-to-Design distances values for different runs in the experimental design. This determines whether the parameter changes in the experimental design result is significant changes in the distance values.
    * signal_validation_oneface.r -- Same as above, without treatment of the face as a variable.
    * viewing_angle_comparison.py -- Directly testing for the influence of the viewing angle on the distance comparisons. **NOT USED IN THE PAPER ANALYSIS**
    * wavelet_normalization.ipynb -- A notebook examining the impact of the unnormalized wavelet tranform, as used in Gill & Henriksen (1990), and the normalized version we, and other wavelet-based transforms like Delta-Variance, use. Deviations from power-law behaviour are more evident in the normalized version.
* The `analysis` folder has functions used for plotting the results and converting the outputs from the cluster jobs.
* `cleanup_scripts` removed old and unused results from the HDF5 data files. There is probably no need to use any of these.
* The `complete` folder has a separate job to pre-compute the dendrograms of the COMPLETE data and save it for future use in order to save computational time. It also has the masking script to produce a signal mask for the observational data. This signal mask was used when creating the moment arrays for these data.
* The `figures` folder has the scripts for reproducing the figures that are not generated from the analysis pipeline. The power-spectrum figure (`ps.py`) requires the raw simulation data, which we have not made available due to storage restrictions. Please contact us if you have interest in the simulation data.
