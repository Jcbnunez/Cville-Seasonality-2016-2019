#!/usr/bin/env bash
#
#SBATCH -J blsmm_env # A single job name for the array
#SBATCH --ntasks-per-node=1 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 8:00:00 ### 20 minutes per job per 10 permutations
#SBATCH --mem 5G
#SBATCH -o /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3

# sbatch --array=1-1000 ~/Overwintering_18_19/temperatureAverage_yearFactor_GLM/aveTemp_yearFactor.lmer.sh
# sacct -j 26666696
# cat /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_26666696_1.err
### SLURM_ARRAY_TASK_ID=1

cd ~/
Rscript Overwintering_18_19/temperatureAverage_yearFactor_GLM/aveTemp_yearFactor.lmer.R ${SLURM_ARRAY_TASK_ID} 1000
