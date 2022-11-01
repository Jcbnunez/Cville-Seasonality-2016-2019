#!/usr/bin/env bash
#
#SBATCH -J run_collect_dat # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00
#SBATCH --mem 2G
#SBATCH -o /scratch/csm6hg/bottleneck/err/out.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/bottleneck/err/out.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load goolf/7.1.0_3.1.4 R/4.0.3
module load gdal geos proj

guide_file=/scratch/csm6hg/bottleneck/model_paramList_fin4_missed

nMax=$( cat $guide_file  | sed '1d' | cut -f2,3 -d',' |  sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d',' )
nMin=$( cat $guide_file  | sed '1d' | cut -f2,3 -d',' |  sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d',' )

echo $nMax $nMin

Rscript \
/scratch/csm6hg/bottleneck/1.Parse.data.for.ABC.R \
${nMax} \
${nMin}

date
echo "done"
