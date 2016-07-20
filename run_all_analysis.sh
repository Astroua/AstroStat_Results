
# Everything!
# Start by copying the output files off of Jasper.

results_dir='/media/eric/Data_3/Astrostat/results/'
scripts_dir='/home/eric/Dropbox/code_development/AstroStat_Results/'

# Sim comparisons
rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/clean_results/*.h5 ${results_dir}/clean/HDF5_files/
rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/noise_same_results/*.h5 ${results_dir}/noisy/HDF5_files/

# Obs comparisons
make_folder_struct=0
if [[ ${make_folder_struct} -eq 1 ]]; then
    mkdir ${results_dir}/noisy/Obs_to_Obs
    mkdir ${results_dir}/noisy/Obs_to_Fid
    mkdir ${results_dir}/noisy/Obs_to_Fid/HDF5
    mkdir ${results_dir}/noisy/Des_to_Obs
    mkdir ${results_dir}/noisy/Des_to_Obs/HDF5
fi

rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/des_to_obs/*.h5 ${results_dir}/noisy/Des_to_Obs/HDF5/
rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/obs_to_fid/*.h5 ${results_dir}/noisy/Obs_to_Fid/HDF5/
rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/obs_to_obs/complete_comparisons.csv ${results_dir}/noisy/Obs_to_Obs/complete_comparisons_new.csv



cd ${results_dir}/clean
# Run the pipeline for the clean and noisy distances
python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/clean/clean mean
# Run the analysis at the common free-fall time.
# python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/clean/clean_freefall freefall

# Noisy
cd ../noisy
python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/noisy/noisy mean

# python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/noisy/noisy_freefall freefall

# Also run the observational analysis from noisy
# Need some way of getting the name of the newly created noise analysis folder (with timestep)
# python ${scripts_dir}/obs_analysis_pipeline.py . NOISE_FOLDER