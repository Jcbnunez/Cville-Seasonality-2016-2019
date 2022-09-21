#!/usr/bin/env bash
#
#SBATCH -J recomp # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 5:00:00 ### 
#SBATCH --mem 40G
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-19


ifile=$(cat files.to.recom.txt | sed "${SLURM_ARRAY_TASK_ID}q;d" )
echo $ifile

#decompress
gzip -d $ifile.gz

#spot check
head $ifile
tail $ifile

#recompress
gzip $ifile


date
echo "done"
