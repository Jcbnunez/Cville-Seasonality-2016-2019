#!/usr/bin/env bash
#
#SBATCH -J fsts.ag # A single job name for the array
#SBATCH --ntasks-per-node=20 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 15:00:00 ### 
#SBATCH --mem 80G
#SBATCH -o slurmOut/fsts.ag%A_%a.out # Standard output
#SBATCH -e slurmOut/fsts.ag%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load gcc/7.1.0
module load openmpi/3.1.4
module load R/4.1.1
module load gdal
module load proj

Rscript \
--vanilla \
3.collect.fst.files.R \
3L

date
echo "done"

