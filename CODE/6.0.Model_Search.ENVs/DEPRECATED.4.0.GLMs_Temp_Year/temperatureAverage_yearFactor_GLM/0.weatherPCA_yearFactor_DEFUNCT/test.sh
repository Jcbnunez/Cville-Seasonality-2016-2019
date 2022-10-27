#!/usr/bin/env bash
#
#SBATCH -J blsmm_env # A single job name for the array
#SBATCH --ntasks-per-node=5 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 0:00:05 ### 20 minutes per job per 10 permutations
#SBATCH --mem 1M
#SBATCH -o /scratch/aob2x/lasso_dest/slurmOut/25869239_%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# sbatch --array=1-100 ~/Overwintering_18_19/temperatureAverage_yearFactor_GLM/test.sh
# sacct -j 25869239
# cat /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_25869239*.out

echo ${SLURM_ARRAY_TASK_ID}
