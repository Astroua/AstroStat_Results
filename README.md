# AstroStat_Results
Scripts for reproducing the results of Koch et al. (2016)

# Running the analysis:

* Jasper scripts -- Main comparisons are run on the Jasper cluster. This includes time-averaged and freefall comparisons of the clean cubes, noise-added cubes, hot fiducials, and observational comparisons. All jobs can be submitted through `jasper/run_all.sh`.
* Smaller analyses -- We do a few more comparisons that are smaller in scope, and so are not setup for running on the cluster. This includes: the effect of AMR (`AMR_comparison_analysis.py`), and the resolution comparison (`resolution_comparison_analysis.py`; we do 2 runs with the 256 grid cubes, and another regridding back down to 128 pixels).
