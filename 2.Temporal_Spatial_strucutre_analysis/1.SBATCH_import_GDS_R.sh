#!/usr/bin/env bash
#
#SBATCH -J MineGDS # A single job name for the array
#SBATCH -c 10
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00 #<= this may depend on your resources
#SBATCH --mem 40G #<= this may depend on your resources
#SBATCH -o ./MineGDS.%A_%a.out # Standard output
#SBATCH -e ./MineGDS.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH -A jcbnunez


module load intel/18.0 intelmpi/18.0
module load goolf/7.1.0_3.1.4
module load gdal proj R/4.0.0

Rscript \
--vanilla \
Import_GDStoR.r \
/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds \
/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv \
DEST.2.0.PoolSNP

