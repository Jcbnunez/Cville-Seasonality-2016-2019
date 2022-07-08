#!/usr/bin/env bash
#
#SBATCH -J run_core20Enrch # A single job name for the array
#SBATCH --ntasks-per-node=4 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 0:60:00 ### 20 minutes per job per 10 permutations
#SBATCH --mem 20G
#SBATCH -o ./slurmOut/out.win.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/out.win.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --array=1-101

module load gcc/7.1.0
module load openmpi/3.1.4
module load R/4.1.1
module load gdal
module load proj

### SLURM_ARRAY_TASK_ID=1
date

Rscript \
./1.0.OmnibusNasaPower.Enrichment.R \
${SLURM_ARRAY_TASK_ID} \
"temp.ave;9;3.Europe_E"

date
