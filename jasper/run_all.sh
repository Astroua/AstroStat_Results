
# Submit jobs to run all comparisons

NODE=1
PROCS=12
PMEM=4000mb
HOURS=72

SCRIPT_PATH=/home/ekoch/code_repos/AstroStat_Results/

# Noiseless
DATA_DIR=/lustre/home/ekoch/sims/SimSuite8/
ADD_NOISE=F
RESULTS_DIR=/lustre/home/ekoch/sims/results/clean_results/

for face1 in {0,2}; do
    for face2 in {0,2}; do
        qsub -N fiducial_comp_"$face1"_"$face2" -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FIDUCIAL="fid_comp",FACE1=$face1,FACE2=$face2,DATA_DIR=$DATA_DIR,ADD_NOISE=$ADD_NOISE,RESULTS_DIR=$RESULTS_DIR $SCRIPT_PATH/jasper/fiducial_submit.pbs
        for fid in {0..4}; do
            qsub -N fiducial_"$fid"_"$face1"_"$face2" -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FIDUCIAL=$fid,FACE1=$face1,FACE2=$face2,DATA_DIR=$DATA_DIR,ADD_NOISE=$ADD_NOISE,RESULTS_DIR=$RESULTS_DIR $SCRIPT_PATH/jasper/fiducial_submit.pbs
        done
    done
done

# Noisy
DATA_DIR=/lustre/home/ekoch/sims/SimSuite8_noise/
ADD_NOISE=T
RESULTS_DIR=/lustre/home/ekoch/sims/results/noise_same_results/

for face1 in {0,2}; do
    for face2 in {0,2}; do
        qsub -N fiducial_noise_comp_"$face1"_"$face2" -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FIDUCIAL="fid_comp",FACE1=$face1,FACE2=$face2,DATA_DIR=$DATA_DIR,ADD_NOISE=$ADD_NOISE,RESULTS_DIR=$RESULTS_DIR $SCRIPT_PATH/jasper/fiducial_submit.pbs
        for fid in {0..4}; do
            qsub -N fiducial_noise_"$fid"_"$face1"_"$face2" -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 -v SCRIPT_PATH=$SCRIPT_PATH,FIDUCIAL=$fid,FACE1=$face1,FACE2=$face2,DATA_DIR=$DATA_DIR,ADD_NOISE=$ADD_NOISE,RESULTS_DIR=$RESULTS_DIR $SCRIPT_PATH/jasper/fiducial_submit.pbs
        done
    done
done

# Obs to Obs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/complete_to_complete.pbs

# Des to Obs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/des_to_complete.pbs

# Obs to Fid
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/complete_to_fid.pbs

