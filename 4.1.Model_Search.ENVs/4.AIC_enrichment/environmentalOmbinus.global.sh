#!/usr/bin/env bash
#
#SBATCH -J blsmm_env # A single job name for the array
#SBATCH --ntasks-per-node=5 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 10:00:00 ### 20 minutes per job per 10 permutations
#SBATCH --mem 40G
#SBATCH -o /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1

# sbatch --array=984-5000%100 /scratch/aob2x/Overwintering_18_19/EnvironmentalOmnibus_Global/environmentalOmbinus.global.sh
# sacct -j 39112815 | tail
# cat /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_38881169_1.err
# head /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_28561442_1.out # Fri Nov 12 10:36:25 EST 2021
# tail /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_28561442_1.out # Fri Nov 12 11:14:06 EST 2021

### SLURM_ARRAY_TASK_ID=1
date

cd /scratch/aob2x/
Rscript Overwintering_18_19/EnvironmentalOmnibus_Global/environmentalOmbinus.global.R ${SLURM_ARRAY_TASK_ID} 5000

date
