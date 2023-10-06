#!/bin/bash
#SBATCH -t 0-00:45:00
#SBATCH -N 1
#SBATCH --mem=15G
#SBATCH -p standard
#SBATCH -c 10
#SBATCH -A berglandlab_standard


module load gcc/9.2.0
module load openmpi/3.1.6
module load R/4.2.1
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p poolcontrol.txt)
Rscript glm_equiv_poolgen.R --args $OPTS
sed -i 's/NA/0/g' *.genobaypass 
