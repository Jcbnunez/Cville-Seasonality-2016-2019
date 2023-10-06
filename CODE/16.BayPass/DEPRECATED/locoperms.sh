#!/bin/bash
#SBATCH -t 4-00:00:00
#SBATCH -N 1
#SBATCH --mem=9G
#SBATCH -p standard
#SBATCH -c 8
#SBATCH -A berglandlab_standard

#WIP

module load gcc/9.2.0
module load openmpi/3.1.6
module load R/4.2.1
Rscript /scratch/nzx3cc/nzx3cc/scripts/permuter.R --args whole_${SLURM_ARRAY_TASK_ID} $SLURM_ARRAY_TASK_ID
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p locostitchcontrol.txt)
baypass="/home/nzx3cc/baypass_2.4/sources/g_baypass"
$baypass -gfile ${OPTS}pool.genobaypass -poolsizefile ${OPTS}pool.poolsize -outprefix "${OPTS}_locostitch" -nthreads 8 -omegafile no${OPTS}_mat_omega.out -efile /scratch/nzx3cc/nzx3cc/perms/whole_${SLURM_ARRAY_TASK_ID}_perm_envtable.txt -contrastfile /scratch/nzx3cc/nzx3cc/perms/whole_${SLURM_ARRAY_TASK_ID}_perm_contrasttable.txt

