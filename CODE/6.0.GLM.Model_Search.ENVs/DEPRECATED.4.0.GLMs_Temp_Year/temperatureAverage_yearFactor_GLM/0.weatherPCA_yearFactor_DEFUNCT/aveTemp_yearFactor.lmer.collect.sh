#!/usr/bin/env bash
#
#SBATCH -J blsmm_env # A single job name for the array
#SBATCH --ntasks-per-node=5 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 0:30:00 ### 20 minutes per job per 10 permutations
#SBATCH --mem 20G
#SBATCH -o /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

module load intel/18.0 intelmpi/18.0 R/3.6.3

# sbatch --array=1-11 ~/Overwintering_18_19/temperatureAverage_yearFactor_GLM/aveTemp_yearFactor.lmer.collect.sh
# sacct -j 26683328
# cat /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_26683189_1.err
### SLURM_ARRAY_TASK_ID=1

cd ~/
Rscript Overwintering_18_19/temperatureAverage_yearFactor_GLM/aveTemp_yearFactor.lmer.collect.R ${SLURM_ARRAY_TASK_ID}
