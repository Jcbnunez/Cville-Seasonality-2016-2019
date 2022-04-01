#!/usr/bin/env bash
#
#SBATCH -J blsmm_env # A single job name for the array
#SBATCH --ntasks-per-node=1 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 0:10:00 ### 10 minutes per job per 10 permutations
#SBATCH --mem 5G
#SBATCH -o /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

module load intel/18.0 intelmpi/18.0 R/3.6.3

# sbatch --array=1-707 ~/Overwintering_18_19/temperatureAverage_yearFactor_GLM/1.3.aveTemp_yearFactor_model_nested_qb/aveTemp_yearFactor.morePerms.collect_part2.sh
# sacct -j 28902106
# cat /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_28888648_1.out
### SLURM_ARRAY_TASK_ID=1

 cd ~/
 Rscript Overwintering_18_19/temperatureAverage_yearFactor_GLM/1.3.aveTemp_yearFactor_model_nested_qb/aveTemp_yearFactor.morePerms.collect.R ${SLURM_ARRAY_TASK_ID}
