#!/usr/bin/env bash
#
#SBATCH -J run_collect_dat # A single job name for the array
#SBATCH -c 2 
#SBATCH -N 1 # on one node
#SBATCH -t 5:00:00 ### 4 hours
#SBATCH --mem 8G
#SBATCH -o ./slurmOut/out.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/out.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-1653

module load gcc/7.1.0
module load openmpi/3.1.4
module load R/4.1.1
module load gdal
module load proj

guide_file=/project/berglandlab/connor/bottleneck/model_paramList_fin3_reduced

nMax=$( cat $guide_file  | sed '1d' | cut -f2,3 -d',' |  sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d',' )
nMin=$( cat $guide_file  | sed '1d' | cut -f2,3 -d',' |  sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d',' )

echo $nMax $nMin  
  
Rscript \
--vanilla \
1.Parse.data.for.ABC.R \
${nMax} \
${nMin} 

date
echo "done"