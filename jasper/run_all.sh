
# Submit jobs to run all comparisons

NODE=1
PROCS=12
PMEM=4000mb
HOURS=72

SCRIPT_PATH=/home/ekoch/code_repos/AstroStat_Results/

# Noiseless
DATA_PATH=/lustre/home/ekoch/sims/SimSuite8/
ADD_NOISE=F
OUTPUT_PATH=/lustre/home/ekoch/sims/results/clean_results/

for face1 in {0,2}; do
    for face2 in {0,2}; do
        qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial_submit.pbs fid_comp $face1 $face2 $DATA_PATH $ADD_NOISE $OUTPUT_PATH
        for fid in {0..4}; do
            qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial_submit.pbs $fid $face1 $face2 $DATA_PATH $ADD_NOISE $OUTPUT_PATH

# Noisy
DATA_PATH=/lustre/home/ekoch/sims/SimSuite8_noise/
ADD_NOISE=T
OUTPUT_PATH=/lustre/home/ekoch/sims/results/noise_same_results/

for face1 in {0,2}; do
    for face2 in {0,2}; do
        qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial_submit.pbs fid_comp $face1 $face2 $DATA_PATH $ADD_NOISE $OUTPUT_PATH
        for fid in {0..4}; do
            qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial_submit.pbs $fid $face1 $face2 $DATA_PATH $ADD_NOISE $OUTPUT_PATH

# Obs to Obs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/complete_to_complete.pbs

# Des to Obs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/des_to_complete.pbs

# Obs to Fid
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/complete_to_fid.pbs

