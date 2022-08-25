#!/usr/bin/env bash
#
#SBATCH -J run_fst_glm # A single job name for the array
#SBATCH -c 2 
#SBATCH -N 1 # on one node
#SBATCH -t 10:00:00 ### 4 hours
#SBATCH --mem 120G
#SBATCH -o ./slurmOut/out.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/out.%A_%a.err # Standard error
#SBATCH -p largemem
#SBATCH --account jcbnunez

module load gcc/7.1.0
module load openmpi/3.1.4
module load R/4.1.1
module load gdal
module load proj

Rscript \
--vanilla \
7.Snp.wise.Analysis.V2.Cville.R

date
echo "done"