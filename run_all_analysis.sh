
# Everything!
# Start by copying the output files off of Jasper.

results_dir='/media/eric/Data_3/Astrostat/results/'
scripts_dir='/home/eric/Dropbox/code_development/AstroStat_Results/'

# Set what to run
make_folder_struct=1
file_copy=1
clean_analysis=1
noisy_analysis=1
obs_analysis=1
ff_analysis=1

# Create the expected folder structure
if [[ ${make_folder_struct} -eq 1 ]]; then
    # Normal analysis
    mkdir ${results_dir}/clean
    mkdir ${results_dir}/clean/HDF5_files
    mkdir ${results_dir}/noisy
    mkdir ${results_dir}/noisy/HDF5_files
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

fi

if [[ ${file_copy} -eq 1 ]]; then
    echo "Copying files from Jasper"
    # Sim comparisons
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/clean_results/*.h5 ${results_dir}/clean/HDF5_files/
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/noise_same_results/*.h5 ${results_dir}/noisy/HDF5_files/

    # Obs comparisons
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/des_to_obs/*.h5 ${results_dir}/noisy/Des_to_Obs/HDF5/
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/obs_to_fid/*.h5 ${results_dir}/noisy/Obs_to_Fid/HDF5/
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/obs_to_obs/complete_comparisons.csv ${results_dir}/noisy/Obs_to_Obs/complete_comparisons_new.csv

    # FF comparisons
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/clean_results_freefall/*.h5 ${results_dir}/clean_freefall/HDF5_files/
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/noise_same_results_freefall/*.h5 ${results_dir}/noisy_freefall/HDF5_files/

fi

# Run the pipeline for the clean and noisy distances

if [[ ${clean_analysis} -eq 1 ]]; then
    echo "Running clean analysis"
    cd ${results_dir}/clean
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/clean/clean mean
fi

# Noisy
if [[ ${noisy_analysis} -eq 1 ]]; then
    echo "Running noisy analysis"
    cd ${results_dir}/noisy
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/noisy/noisy mean
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
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/clean_freefall/clean_freefall mean
    echo "Running noisy freefall analysis"
    cd ${results_dir}/noisy_freefall
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/noisy_freefall/noisy/noisy_freefall mean
fi

