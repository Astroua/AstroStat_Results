
# Submit jobs to run all comparisons

NODE=1
PROCS=12
PMEM=4gb
HOURS=72

SCRIPT_PATH=/home/ekoch/code_repos/TurbuStat/Examples/

# Noiseless
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial0_all.pbs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial1_all.pbs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial2_all.pbs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial3_all.pbs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial4_all.pbs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial_comp.pbs

# Noisy
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial0_all_noise.pbs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial1_all_noise.pbs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial2_all_noise.pbs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial3_all_noise.pbs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial4_all_noise.pbs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/fiducial_comp_noise.pbs

# Obs to Obs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/complete_to_complete.pbs

# Des to Obs
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/des_to_complete.pbs

# Obs to Fid
qsub -l nodes=$NODE:ppn=$PROCS,pmem=$PMEM,walltime=$HOURS:00:00 $SCRIPT_PATH/jasper/complete_to_fid.pbs

