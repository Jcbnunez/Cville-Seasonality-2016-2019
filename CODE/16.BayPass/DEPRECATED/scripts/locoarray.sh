#!/bin/bash
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem 9G
#SBATCH -t 6-20:30:00
#SBATCH -p standard
#SBATCH -A berglandlab_standard

module load gcc/9.2.0
module load openmpi/3.1.6
module load R/4.2.1
baypass="/home/nzx3cc/baypass_2.4/sources/g_baypass"
  Rscript "/scratch/nzx3cc/nzx3cc/scripts/poolextractionLOCO_${SLURM_ARRAY_TASK_ID}.R" --args 2000 "lococomparison_${SLURM_ARRAY_TASK_ID}"
  sed -i 's/NA/0/g' lococomparison_${SLURM_ARRAY_TASK_ID}.genobaypass
  $baypass -gfile lococomparison_${SLURM_ARRAY_TASK_ID}.genobaypass -poolsizefile lococomparison_${SLURM_ARRAY_TASK_ID}.poolsize -outprefix "lococomparison_${SLURM_ARRAY_TASK_ID}" -nthreads 8 

 

  
