#!/usr/bin/env bash
#
#SBATCH -J aic.test # A single job name for the array
#SBATCH --ntasks-per-node=20 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 
#SBATCH --mem 80G
#SBATCH -o slurmOut/test.ag%A_%a.out # Standard output
#SBATCH -e slurmOut/test.ag%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-2678

module load gcc/7.1.0
module load openmpi/3.1.4
module load R/4.1.1
module load gdal
module load proj

Rscript \
--vanilla \
1.0.PLOT.AIC.delta.suppFig.R \
${SLURM_ARRAY_TASK_ID}

date
echo "done"

