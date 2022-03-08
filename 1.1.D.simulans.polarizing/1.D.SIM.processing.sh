#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=8:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez

### Include D. Simulans data for polarizing

module load bcftools
module load tabix

###
bcftools view -s "SIM" /project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.vcf.gz > \
SIM.dmel6.vcf

bgzip SIM.dmel6.vcf
tabix SIM.dmel6.vcf.gz


