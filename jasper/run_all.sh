
# Submit jobs to run all comparisons

NODE=1
PROCS=12
PMEM=4000mb
HOURS=72

# Set which analyses to run
run_clean=1
run_noise=1
run_regrid=1
run_obs=1
run_fid_reproc=1
run_face_compare=1
run_hotfid_compare=1

SCRIPT_PATH=/home/ekoch/code_repos/AstroStat_Results/

# Set the type of comparison to run. If 'max', all timesteps are run against
# the timesteps of all others (starting from the last timestep in each).
# If 'freefall', comparisons are only done between the timesteps in each simulation
# nearest to a free-fall time.
# COMPARE_TYPE='max'
# COMPARE_TYPE=freefall

for COMPARE_TYPE in 'max', 'freefall'; do

    echo "Running comparison type: "$COMPARE_TYPE

    # Noiseless
    DATA_DIR=/lustre/home/ekoch/sims/SimSuite8/
    ADD_NOISE=F
    if [[ $COMPARE_TYPE = max ]]; then
        RESULTS_DIR=/lustre/home/ekoch/sims/results/clean_results/
    else
        RESULTS_DIR=/lustre/home/ekoch/sims/results/clean_results_freefall/
    fi

    if [[ $run_clean -eq 1 ]]; then
        # Run face 0 and face 2 design to fiducial comparisons
        for face1 in {0,2}; do
            for face2 in {0,2}; do
                for fid in {0..4}; do
                    qsub -N fiducial_"$fid"_"$face1"_"$face2"_"$COMPARE_TYPE" -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FIDUCIAL=$fid,FACE_1=$face1,FACE_2=$face2,DATA_DIR=$DATA_DIR,ADD_NOISE=$ADD_NOISE,RESULTS_DIR=$RESULTS_DIR,COMPARE_TYPE=$COMPARE_TYPE $SCRIPT_PATH/jasper/fiducial_submit.pbs
                done
            done
        done

        # Run all faces of the fiducial comparisons
        for face1 in {0..2}; do
            for face2 in {0..2}; do
                qsub -N fiducial_comp_"$face1"_"$face2"_"$COMPARE_TYPE" -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FIDUCIAL="fid_comp",FACE_1=$face1,FACE_2=$face2,DATA_DIR=$DATA_DIR,ADD_NOISE=$ADD_NOISE,RESULTS_DIR=$RESULTS_DIR,COMPARE_TYPE=$COMPARE_TYPE $SCRIPT_PATH/jasper/fiducial_submit.pbs
            done
        done
    fi

    # Noisy
    DATA_DIR=/lustre/home/ekoch/sims/SimSuite8_noise/
    ADD_NOISE=T
    if [[ $COMPARE_TYPE = max ]]; then
        RESULTS_DIR=/lustre/home/ekoch/sims/results/noise_same_results/
    else
        RESULTS_DIR=/lustre/home/ekoch/sims/results/noise_same_results_freefall/
    fi

    if [[ $run_noise -eq 1 ]]; then
        for face1 in {0,2}; do
            for face2 in {0,2}; do
                qsub -N fiducial_noise_comp_"$face1"_"$face2"_"$COMPARE_TYPE" -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FIDUCIAL="fid_comp",FACE_1=$face1,FACE_2=$face2,DATA_DIR=$DATA_DIR,ADD_NOISE=$ADD_NOISE,RESULTS_DIR=$RESULTS_DIR,COMPARE_TYPE=$COMPARE_TYPE $SCRIPT_PATH/jasper/fiducial_submit.pbs
                for fid in {0..4}; do
                    qsub -N fiducial_noise_"$fid"_"$face1"_"$face2"_"$COMPARE_TYPE" -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FIDUCIAL=$fid,FACE_1=$face1,FACE_2=$face2,DATA_DIR=$DATA_DIR,ADD_NOISE=$ADD_NOISE,RESULTS_DIR=$RESULTS_DIR,COMPARE_TYPE=$COMPARE_TYPE $SCRIPT_PATH/jasper/fiducial_submit.pbs
                done
            done
        done
    fi

    # Regridded
    DATA_DIR=/lustre/home/ekoch/sims/SimSuite8_regrid/
    ADD_NOISE=T
    if [[ $COMPARE_TYPE = max ]]; then
        RESULTS_DIR=/lustre/home/ekoch/sims/results/regrid_results/
    else
        RESULTS_DIR=/lustre/home/ekoch/sims/results/regrid_results_freefall/
    fi

    if [[ $run_regrid -eq 1 ]]; then
        for face1 in {0,2}; do
            for face2 in {0,2}; do
                qsub -N fiducial_regrid_comp_"$face1"_"$face2"_"$COMPARE_TYPE" -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FIDUCIAL="fid_comp",FACE_1=$face1,FACE_2=$face2,DATA_DIR=$DATA_DIR,ADD_NOISE=$ADD_NOISE,RESULTS_DIR=$RESULTS_DIR,COMPARE_TYPE=$COMPARE_TYPE $SCRIPT_PATH/jasper/fiducial_submit.pbs
                for fid in {0..4}; do
                    qsub -N fiducial_regrid_"$fid"_"$face1"_"$face2"_"$COMPARE_TYPE" -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FIDUCIAL=$fid,FACE_1=$face1,FACE_2=$face2,DATA_DIR=$DATA_DIR,ADD_NOISE=$ADD_NOISE,RESULTS_DIR=$RESULTS_DIR,COMPARE_TYPE=$COMPARE_TYPE $SCRIPT_PATH/jasper/fiducial_submit.pbs
                done
            done
        done
    fi
done

# Only run the observational comparisons when using 'max' (to avoid re-computing)
if [[ $run_obs -eq 1 ]]; then
    # Obs to Obs
    qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH $SCRIPT_PATH/jasper/complete_to_complete.pbs

    # Des to Obs
    qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH $SCRIPT_PATH/jasper/des_to_complete.pbs

    # Obs to Fid
    qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH $SCRIPT_PATH/jasper/complete_to_fid.pbs

fi


# Run reprocessed fiducial comparison (for resolution analysis)
DATA_DIR=/lustre/home/ekoch/sims/Fiducial_reproc/
ADD_NOISE=F
RESULTS_DIR=/lustre/home/ekoch/sims/results/fiducial_reproc

if [[ $run_fid_reproc -eq 1 ]]; then
    face1=0
    face2=0
    qsub -N reproc_fiducial_noise_comp_"$face1"_"$face2" -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FIDUCIAL="fid_comp",FACE_1=$face1,FACE_2=$face2,DATA_DIR=$DATA_DIR,ADD_NOISE=$ADD_NOISE,RESULTS_DIR=$RESULTS_DIR,COMPARE_TYPE=$COMPARE_TYPE $SCRIPT_PATH/jasper/fiducial_submit.pbs
fi

DATA_DIR=/lustre/home/ekoch/sims/SimSuite8/
RESULTS_DIR=/lustre/home/ekoch/sims/results/face_compare/
if [[ $run_face_compare -eq 1 ]]; then
    # Face 0 to Face 1
    qsub -N face_comparison_0_1 -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FACE_1=0,FACE_2=1,DATA_DIR=$DATA_DIR,RESULTS_DIR=$RESULTS_DIR $SCRIPT_PATH/jasper/run_viewing_angle_distances.pbs
    # Face 0 to Face 2
    qsub -N face_comparison_0_2 -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FACE_1=0,FACE_2=2,DATA_DIR=$DATA_DIR,RESULTS_DIR=$RESULTS_DIR $SCRIPT_PATH/jasper/run_viewing_angle_distances.pbs
    # Face 1 to Face 2
    qsub -N face_comparison_1_2 -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FACE_1=1,FACE_2=2,DATA_DIR=$DATA_DIR,RESULTS_DIR=$RESULTS_DIR $SCRIPT_PATH/jasper/run_viewing_angle_distances.pbs
fi

DATA_DIR=/lustre/home/ekoch/sims/SimSuite8_hot/
RESULTS_DIR=/lustre/home/ekoch/sims/results/hot_compare/
if [[ $run_hotfid_compare -eq 1 ]]; then
    # Run Fiducial comparisons to the hot Fiducials
    # The clean normal fiducial to fiducial comparisons are valid here, so no need to re-run
    for face1 in {0,2}; do
        for face2 in {0,2}; do
            # qsub -N fiducial_regrid_comp_"$face1"_"$face2" -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FIDUCIAL="fid_comp",FACE_1=$face1,FACE_2=$face2,DATA_DIR=$DATA_DIR,ADD_NOISE=$ADD_NOISE,RESULTS_DIR=$RESULTS_DIR,COMPARE_TYPE=$COMPARE_TYPE $SCRIPT_PATH/jasper/fiducial_submit.pbs
            for fid in {0..4}; do
                qsub -N fiducial_hot_"$fid"_"$face1"_"$face2" -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FIDUCIAL=$fid,FACE_1=$face1,FACE_2=$face2,DATA_DIR=$DATA_DIR,ADD_NOISE=$ADD_NOISE,RESULTS_DIR=$RESULTS_DIR,COMPARE_TYPE=$COMPARE_TYPE,HOT_RUN=T $SCRIPT_PATH/jasper/fiducial_submit.pbs
            done
        done
    done
fi
