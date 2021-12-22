#!/usr/bin/env bash
#
#SBATCH -J blsmm_env # A single job name for the array
#SBATCH --ntasks-per-node=1 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 2:00:00 ### 20 minutes per job per 10 permutations
#SBATCH --mem 5G
#SBATCH -o /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


module load intel/18.0 intelmpi/18.0 R/3.6.3

# sbatch --array=1-1000 ~/Overwintering_18_19/weatherPCA_yearFactor/aveTemp_yearFactor.sh
# sacct -j 25647087
# cat /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_25647087_500.err
### SLURM_ARRAY_TASK_ID=1

cd ~/
Rscript Overwintering_18_19/weatherPCA_yearFactor/aveTemp_yearFactor.R ${SLURM_ARRAY_TASK_ID} 1000
