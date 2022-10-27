#!/usr/bin/env bash
#
#SBATCH -J blsmm_env # A single job name for the array
#SBATCH --ntasks-per-node=1 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 20 minutes per job per 10 permutations
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

module load parallel intel/18.0 intelmpi/18.0 R/3.6.3

# sbatch --array=1-7 ~/Overwintering_18_19/temperatureAverage_yearFactor_GLM/aveTemp_yearFactor.morePerms.collect_part1.sh
# sacct -j 28719915
# cat /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_28630032_5.err
### SLURM_ARRAY_TASK_ID=1

# cd ~/
# Rscript Overwintering_18_19/temperatureAverage_yearFactor_GLM/aveTemp_yearFactor.morePerms.collect.R ${SLURM_ARRAY_TASK_ID}

# SLURM_ARRAY_TASK_ID=1

dir=$( ls -d /scratch/aob2x/dest_glm_morePerms/* | sed "${SLURM_ARRAY_TASK_ID}q;d" )
echo $dir

split_fun () {
  fn=${1} # fn="/scratch/aob2x/dest_glm_morePerms/DE_Bro/DE_Bro;allPerm;aveTemp_yearFactor;100.csv"
  echo ${fn}
  cd /scratch/aob2x/dest_morePerms
  cat ${fn} |
  awk -F ',' '{
    print $0 >> $2"_"$3".csv"
  }'
}
export -f split_fun

parallel -j1 split_fun ::: $( ls -d $dir/* )
