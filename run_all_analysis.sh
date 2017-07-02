
# Everything!
# Start by copying the output files off of Jasper.

results_dir='/home/eric/Dropbox/Data/Astrostat'
scripts_dir='/home/eric/Dropbox/code_development/AstroStat_Results/'
data_path='/mnt/MyRAID/Astrostat'

# Which face(s) to analyze
faces="0"

# Set what to run
make_folder_struct=1
file_copy=1
clean_analysis=1
noisy_analysis=1
regrid_analysis=1
obs_analysis=1
ff_analysis=1
res_comp_analysis=1
face_comparison=1
hot_comparison=1
basic_property_comparison=1
amr_comparison=1

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
    # mkdir /media/eric/MyRAID/Astrostat/Fiducial_256/moments
    # mkdir /media/eric/MyRAID/Astrostat/Fiducial_reproc/moments
    # Face comparison
    mkdir ${results_dir}/face_compare
    # Hot fiducial comparison
    mkdir ${results_dir}/hot_compare
    mkdir ${results_dir}/hot_compare/HDF5_files

    # Basic property comparison between cubes
    mkdir ${results_dir}/basic_property_comparison

    # Add the design file to the results dir
    cp ~/Dropbox/AstroStatistics/Full_Factorial/Design7Matrix.csv ${results_dir}

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
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/obs_to_obs/complete_comparisons.csv ${results_dir}/noisy/Obs_to_Obs/

    # FF comparisons
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/clean_results_freefall/*.h5 ${results_dir}/clean_freefall/HDF5_files/
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/noise_same_results_freefall/*.h5 ${results_dir}/noisy_freefall/HDF5_files/
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/regrid_results_freefall/*.h5 ${results_dir}/regrid_freefall/HDF5_files/

    # Reproc fiducial comparisons
    # rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/fiducial_reproc/*.h5 ${results_dir}/resolution_comparison

    # Face comparisons
    # rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/face_compare/*.csv ${results_dir}/face_compare

    # Hot fiducial comparisons
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/hot_compare/*.h5 ${results_dir}/hot_compare/HDF5_files/
    # Need the fiducial-fiducial clean comparisons here (no need to be re-running them)
    rsync ekoch@jasper.westgrid.ca:/home/ekoch/sims/results/clean_results/SimSuite8_fiducialfid_comp*.h5 ${results_dir}/hot_compare/HDF5_files/

fi

# Run the pipeline for the clean and noisy distances

if [[ ${clean_analysis} -eq 1 ]]; then
    echo "Running clean analysis"
    cd ${results_dir}/clean
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/clean/clean mean F ${faces}
fi

# Noisy
if [[ ${noisy_analysis} -eq 1 ]]; then
    echo "Running noisy analysis"
    cd ${results_dir}/noisy
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/noisy/noisy mean T ${faces}
fi

# Regrid
if [[ ${regrid_analysis} -eq 1 ]]; then
    echo "Running regrid analysis"
    cd ${results_dir}/regrid
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/regrid/regrid mean F ${faces}
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
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/clean_freefall/clean_freefall mean F ${faces}
    echo "Running noisy freefall analysis"
    cd ${results_dir}/noisy_freefall
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/noisy_freefall/noisy_freefall mean T ${faces}
    echo "Running regrid freefall analysis"
    cd ${results_dir}/regrid_freefall
    python ${scripts_dir}/analysis_pipeline.py . ${results_dir}/Design7Matrix.csv ${scripts_dir} ${results_dir}/regrid_freefall/regrid_freefall mean F ${faces}
fi

if [[ ${res_comp_analysis} -eq 1 ]]; then
    echo "Running resolution comparison"
    create_moments=0
    # Need to create the moments first.
    if [[ ${create_moments} -eq 1 ]]; then
        python ${scripts_dir}/jasper/reduce_and_save_moments.py ${data_path}/Fiducial_256/ ${data_path}/Fiducial_256/moments/ T F
        python ${scripts_dir}/jasper/reduce_and_save_moments.py ${data_path}/Fiducial_reproc/ ${data_path}/Fiducial_128_reproc/moments/ T F
    fi

    # Run on w/ the full 256
    python ${scripts_dir}/resolution_comparison_analysis.py F ${results_dir}
    # And with the 256 regridded to 128
    python ${scripts_dir}/resolution_comparison_analysis.py T ${results_dir}
fi

if [[ ${face_comparison} -eq 1 ]]; then
    python ${scripts_dir}/viewing_angle_comparison.py ${results_dir}/face_compare
fi

if [[ ${hot_comparison} -eq 1 ]]; then
    cd ${results_dir}/hot_compare
    python ${scripts_dir}/hot_compare_analysis.py . ${scripts_dir} ${results_dir}/hot_compare/hot_compare mean F ${faces}
fi

if [[ ${basic_property_comparison} -eq 1 ]]; then

    # Compare peak flux, total flux, and average linewidth properties
    # Original sim set
    echo "Running basic property comparisons."
    python ${scripts_dir}/model_basic_properties.py ${data_path}/SimSuite8/ ${results_dir}/Design7Matrix.csv ${results_dir}/basic_property_comparison
    # Regridded sim set
    python ${scripts_dir}/model_basic_properties.py ${data_path}/SimSuite8_regrid/ ${results_dir}/Design7Matrix.csv ${results_dir}/basic_property_comparison

fi

if [[ ${amr_comparison} -eq 1 ]]; then

    echo "Running AMR comparison"
    python ${scripts_dir}/AMR_comparison_analysis ${data_path}/SimSuite8/ ${data_path}/Fiducial_AMR/ ${results_dir}
fi