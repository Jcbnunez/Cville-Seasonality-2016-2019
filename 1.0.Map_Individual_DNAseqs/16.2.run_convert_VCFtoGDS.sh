#!/usr/bin/env bash
#
#SBATCH -J annotate # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00 ### 6 hours
#SBATCH --mem 20G
#SBATCH -p standard
#SBATCH --account jcbnunez

module load htslib bcftools
module load intel/18.0 intelmpi/18.0
module load goolf/7.1.0_3.1.4
module load gdal proj R/4.0.0
module load tabix

##
echo "make GDS"
Rscript --vanilla ./16.1.convert.vcftogds.r \
CM_pops.AllChrs.Whatshap.shapeit.vcf