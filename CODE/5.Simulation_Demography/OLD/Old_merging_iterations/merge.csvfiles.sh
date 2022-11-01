#!/usr/bin/env bash
#SBATCH -J mega_merge # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH --cpus-per-task=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 #  hours
#SBATCH --mem 5G
#SBATCH -o /scratch/csm6hg/bottleneck/err/bottle.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/bottleneck/err/bottle.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

SLURM_ARRAY_TASK_ID=1
rep=1

# Working directory
wd="/scratch/csm6hg/bottleneck"

# Start job
echo "Start"
date

# Parameter file
paramFile=${wd}/model_paramList_fin3

# Extract constants from parameter file
slurmID=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f1 )
nMax=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
nMin=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )
Rep=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )
nSamp=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f5 )
Gen=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f6 )

# Output directory
output1="/project/berglandlab/connor/bottleneck/"
output="/scratch/csm6hg/bottleneck/tmp"

# Check size of files
SIZE=$(ls -ltr ${wd}/data.condense/freq.output.new.pca.${nMax}.${nMin}.csv \
| tr -s ' ' | cut -d ' ' -f 5)
ZERO=0

# For loop replicate
for rep in {1..100}; do

  # Start replicate
  echo $rep

  # Dependent on successful run
  if [ $SIZE -eq $ZERO ]
  then

    echo $SIZE
    echo "File is equal to zero"

    # Send to scratch
    find ${output1}/data3.pca/ -name "*${rep}.${nMax}.${nMin}*" -exec cp {} ${output} \;
    find ${output1}/data3.afvar/ -name "*${rep}.${nMax}.${nMin}*" -exec cp {} ${output} \;
    find ${output1}/data3.euc/ -name "*${rep}.${nMax}.${nMin}*" -exec cp {} ${output} \;
    find ${output1}/data3.fst/ -name "*${rep}.${nMax}.${nMin}*" -exec cp {} ${output} \;

    # Remove header and rbind PCA
    awk 'FNR > 1' ${output1}/data3.pca/freq.output.new.pca.${rep}.${nMax}.${nMin}.*.csv >> \
    ${wd}/data.reps/freq.output.new.pca.${rep}.${nMax}.${nMin}.csv

    # Remove header and rbind AFVAR
    awk 'FNR > 1' ${output}/freq.output.new.afvar.${nMax}.${nMin}.*.csv > \
    ${wd}/data.condense/freq.output.new.afvar.${nMax}.${nMin}.csv

    # Remove header and rbind EUC
    awk 'FNR > 1' ${output}/freq.output.new.euc.${nMax}.${nMin}.*.csv > \
    ${wd}/data.condense/freq.output.new.euc.${nMax}.${nMin}.csv

    # Remove header and rbind FST
    awk 'FNR > 1' ${output}/freq.output.new.fst.${nMax}.${nMin}.*.csv > \
    ${wd}/data.condense/freq.output.new.fst.${nMax}.${nMin}.csv

    # Remove replicate files
    rm ${output}/freq.output.new.pca.${nMax}.${nMin}.*.csv
    rm ${output}/freq.output.new.afvar.${nMax}.${nMin}.*.csv
    rm ${output}/freq.output.new.euc.${nMax}.${nMin}.*.csv
    rm ${output}/freq.output.new.fst.${nMax}.${nMin}.*.csv

  # Finish
  fi

# Finish
echo "Finish" ${nMax} ${nMin}
date
done

# Check if all files are empty
# find *pca* -type f -size 0 | wc -l
