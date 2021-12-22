#! /bin/bash

#SBATCH -J correlation_analysis # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 3:00:00 #<= this may depend on your resources
#SBATCH --mem 4G #<= this may depend on your resources
#SBATCH -o ./slurmOut/correlation_analysis.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/correlation_analysis.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH -A jcbnunez
#SBATCH --array=1-8713

## 1-8713

module load intel/18.0 intelmpi/18.0
module load goolf/7.1.0_3.1.4
module load gdal proj R/4.0.0

mkdir corrOut

Rscript \
--vanilla \
2.run_correlation_iteration.R \
${SLURM_ARRAY_TASK_ID}

