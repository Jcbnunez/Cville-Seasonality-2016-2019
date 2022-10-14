#!/usr/bin/env bash
#
#SBATCH -J remake.BPGP # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 5:00:00 ### 
#SBATCH --mem 40G
#SBATCH -o slurmOut/glm.fsts.ag%A_%a.out # Standard output
#SBATCH -e slurmOut/glm.fsts.ag%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-201

#find "full.vcf.bz2" ./ -print | grep "full.vcf.bz2"  > files.to.merge.list

i=$( cat files.to.merge.list| sed "${SLURM_ARRAY_TASK_ID}q;d" )
echo $i 

bzip2 -d $i 

echo "done"
date