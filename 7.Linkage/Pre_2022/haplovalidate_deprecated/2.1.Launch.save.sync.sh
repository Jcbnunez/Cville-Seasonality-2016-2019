#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --partition=largemem
#SBATCH --account=jcbnunez

module load intel/18.0 intelmpi/18.0
module load goolf/7.1.0_3.1.4
module load gdal proj R/4.0.0

Rscript \
--vanilla \
2.save_sync_file.r