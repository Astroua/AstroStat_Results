
# Everything!
# Start by copying the output files off of Jasper.

results_dir='/media/eric/Data_3/Astrostat/results/'
scripts_dir='/home/eric/Dropbox/code_development/AstroStat_Results/'

# Set what to run
make_folder_struct=1
file_copy=1
clean_analysis=1
noisy_analysis=1
regrid_analysis=1
obs_analysis=1
ff_analysis=1
res_comp_analysis=1

# Create the expected folder structure
if [[ ${make_folder_struct} -eq 1 ]]; then
    # Normal analysis
    mkdir ${results_dir}/clean
    mkdir ${results_dir}/clean/HDF5_files
    mkdir ${results_dir}/noisy
    mkdir ${results_dir}/noisy/HDF5_files
    mkdir ${results_dir}/regrid
    mkdir ${results_dir}/regrid/HDF5_files
    # Observational comparison analysis
    mkdir ${results_dir}/noisy/Obs_to_Obs
    mkdir ${results_dir}/noisy/Obs_to_Fid
    mkdir ${results_dir}/noisy/Obs_to_Fid/HDF5
    mkdir ${results_dir}/noisy/Des_to_Obs
    mkdir ${results_dir}/noisy/Des_to_Obs/HDF5
    # Common FF analysis
    mkdir ${results_dir}/clean_freefall
    mkdir ${results_dir}/clean_freefall/HDF5_files
    mkdir ${results_dir}/noisy_freefall
    mkdir ${results_dir}/noisy_freefall/HDF5_files
    mkdir ${results_dir}/regrid_freefall
    mkdir ${results_dir}/regrid_freefall/HDF5_files
    # Res comparison
    mkdir /media/eric/Data_3/Astrostat/Fiducial_256/moments
    mkdir /media/eric/Data_3/Astrostat/Fiducial_reproc/moments

fi

if [[ ${file_copy} -eq 1 ]]; then
    echo "Copying files from Jasper"
    # Sim comparisons
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/clean_results/*.h5 ${results_dir}/clean/HDF5_files/
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/noise_same_results/*.h5 ${results_dir}/noisy/HDF5_files/
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/regrid_results/*.h5 ${results_dir}/regrid/HDF5_files/

    # Obs comparisons
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/des_to_obs/*.h5 ${results_dir}/noisy/Des_to_Obs/HDF5/
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/obs_to_fid/*.h5 ${results_dir}/noisy/Obs_to_Fid/HDF5/
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/obs_to_obs/complete_comparisons.csv ${results_dir}/noisy/Obs_to_Obs/complete_comparisons_new.csv

    # FF comparisons
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/clean_results_freefall/*.h5 ${results_dir}/clean_freefall/HDF5_files/
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/noise_same_results_freefall/*.h5 ${results_dir}/noisy_freefall/HDF5_files/
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/regrid_results_freefall/*.h5 ${results_dir}/regrid_freefall/HDF5_files/

    # Reproc fiducial comparisons
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/fiducial_reproc/*.h5 ${results_dir}/resolution_comparison

fi

# Run the pipeline for the clean and noisy distances

if [[ ${clean_analysis} -eq 1 ]]; then
    echo "Running clean analysis"
    cd ${results_dir}/clean
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/clean/clean mean F
fi

# Noisy
if [[ ${noisy_analysis} -eq 1 ]]; then
    echo "Running noisy analysis"
    cd ${results_dir}/noisy
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/noisy/noisy mean T
fi

# Regrid
if [[ ${regrid_analysis} -eq 1 ]]; then
    echo "Running regrid analysis"
    cd ${results_dir}/regrid
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/regrid/regrid mean F
fi

# Also run the observational analysis from noisy
# Need some way of getting the name of the newly created noise analysis folder (with timestep)
# if [[ ${obs_analysis} -eq 0 ]]; then
#     echo "Running observational analysis"
#     python ${scripts_dir}/obs_analysis_pipeline.py . NOISE_FOLDER
# fi

# Run the analysis at the common free-fall time.
if [[ ${ff_analysis} -eq 1 ]]; then
    echo "Running clean freefall analysis"
    cd ${results_dir}/clean_freefall
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/clean_freefall/clean_freefall mean F
    echo "Running noisy freefall analysis"
    cd ${results_dir}/noisy_freefall
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/noisy_freefall/noisy_freefall mean T
    echo "Running regrid freefall analysis"
    cd ${results_dir}/regrid_freefall
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/regrid_freefall/regrid_freefall mean F
fi

if [[ ${res_comp_analysis} -eq 1 ]]; then
    echo "Running resolution comparison"
    create_moments=0
    # Need to create the moments first.
    if [[ ${create_moments} -eq 1 ]]; then
        python ${scripts_dir}/jasper/reduce_and_save_moments.py /media/eric/Data_3/Astrostat/Fiducial_256/ /media/eric/Data_3/Astrostat/Fiducial_256/moments/ T F
        python ${scripts_dir}/jasper/reduce_and_save_moments.py /media/eric/Data_3/Astrostat/Fiducial_reproc/ /media/eric/Data_3/Astrostat/Fiducial_128_reproc/moments/ T F
    fi

    # Run on w/ the full 256
    python ${scripts_dir}/resolution_comparison_analysis.py F
    # And with the 256 regridded to 128
    python ${scripts_dir}/resolution_comparison_analysis.py T
fi

