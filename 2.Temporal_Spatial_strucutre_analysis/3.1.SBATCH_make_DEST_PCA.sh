#!/usr/bin/env bash
#
#SBATCH -J MineGDS # A single job name for the array
#SBATCH -c 2
#SBATCH -N 1 # on one node
#SBATCH -t 10:00:00 #<= this may depend on your resources
#SBATCH --mem 120G #<= this may depend on your resources
#SBATCH -o ./MakePCA.%A_%a.out # Standard output
#SBATCH -e ./MakePCA.%A_%a.err # Standard error
#SBATCH -p largemem
#SBATCH -A jcbnunez


module load intel/18.0 intelmpi/18.0
module load goolf/7.1.0_3.1.4
module load gdal proj R/4.0.0

Rscript \
--vanilla \
4.Make_PCA_time_space_DEST.R



