#!/usr/bin/env bash
#
#SBATCH -J blsmm_env # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 0:120:00 ### 20 minutes per job per 10 permutations
#SBATCH --mem 50G
#SBATCH -o /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1

# sbatch --array=1-101 ~/Overwintering_18_19/pheno/pheno.sh
# sacct -j 32323504
# cat /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_32323504_1.err
# head /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_28561442_1.out # Fri Nov 12 10:36:25 EST 2021
# tail /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_28561442_1.out # Fri Nov 12 11:14:06 EST 2021

### SLURM_ARRAY_TASK_ID=1
date

cd ~/
Rscript Overwintering_18_19/pheno/pheno.R ${SLURM_ARRAY_TASK_ID}

date
