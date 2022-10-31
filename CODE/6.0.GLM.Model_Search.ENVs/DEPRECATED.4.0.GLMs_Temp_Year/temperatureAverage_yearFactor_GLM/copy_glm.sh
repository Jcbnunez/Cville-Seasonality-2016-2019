#!/usr/bin/env bash
#
#SBATCH -J blsmm_env # A single job name for the array
#SBATCH --ntasks-per-node=1 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 4:00:00 ### 20 minutes per job per 10 permutations
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

### sbatch ~/Overwintering_18_19/temperatureAverage_yearFactor_GLM/copy_glm.sh
## sacct -j 29719504
### cat /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_29719504.

# mkdir /project/berglandlab/thermal_glm_dest

rsync \
--recursive \
/scratch/aob2x/dest_glm_morePerms_nested_qb/processedGLM \
/project/berglandlab/thermal_glm_dest/dest_glm_final_nested_qb/.

rsync \
--recursive \
/scratch/aob2x/dest_glm_morePerms_nested_qb/windowAnalysis \
/project/berglandlab/thermal_glm_dest/dest_glm_final_nested_qb/.
