#!/usr/bin/env bash
#
#SBATCH -J recomp # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 8:00:00 ### 
#SBATCH --mem 40G
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-20

ifile=$(cat ./3.files.to.recom.txt | sed "${SLURM_ARRAY_TASK_ID}q;d" )
echo $ifile

### order matters
cd folder_${SLURM_ARRAY_TASK_ID}

#spot check
echo "al head"
head $ifile.al
echo "al tail"
tail $ifile.al
echo "un head"
head $ifile.un
echo "un tail"
tail $ifile.un

### merge --> 
cat $ifile.al $ifile.un > $ifile

#recompress
gzip $ifile

date
echo "done"
