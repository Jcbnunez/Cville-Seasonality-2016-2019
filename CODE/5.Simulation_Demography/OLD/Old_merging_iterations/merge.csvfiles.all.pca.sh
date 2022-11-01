#!/usr/bin/env bash
#SBATCH -J mega_merge # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH --cpus-per-task=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-05:00 #  hours
#SBATCH --mem 1G
#SBATCH -o /scratch/csm6hg/bottleneck/err/bottle.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/bottleneck/err/bottle.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SLURM_ARRAY_TASK_ID=738

# Working directory
wd="/scratch/csm6hg/bottleneck"

# Parameter file
paramFile=${wd}/model_paramList_fin3_reduced

# Extract constants from parameter file
slurmID=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f1 )
nMax=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
nMin=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )
Rep=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )
nSamp=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f5 )
Gen=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f6 )

# Output directory
output="/project/berglandlab/connor/bottleneck/data3"

# Remove header and rbind PCA
awk 'FNR > 1' ${output}/freq.output.new.pca.${nMax}.${nMin}.*.csv > \
/project/berglandlab/connor/bottleneck/data.condense/freq.output.new.pca.${nMax}.${nMin}.csv

# Remove replicate files
rm ${output}/freq.output.new.pca.${nMax}.${nMin}.*.csv

# Finish
echo "Finish" ${nMax} ${nMin}
echo date
