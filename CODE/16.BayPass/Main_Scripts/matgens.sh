#!/bin/bash
#SBATCH -t 4-00:00:00
#SBATCH -A berglandlab_standard
#SBATCH -c 8
#SBATCH -p standard
#SBATCH -N 1
#SBATCH --mem 5G 

module load gcc/9.2.0
module load openmpi/3.1.6
module load R/4.2.1
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p matcontrol.txt)
baypass="/home/nzx3cc/baypass_2.4/sources/g_baypass"
$baypass -gfile ${OPTS}.genobaypass -poolsizefile ${OPTS}.poolsize -outprefix $OPTS -nthreads 8
