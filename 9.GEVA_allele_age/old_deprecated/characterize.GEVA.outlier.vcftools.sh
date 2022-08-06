#### ---- extract genotypes
#### 
#### 

#! /bin/bash

#SBATCH -J hap.count 
#SBATCH --ntasks-per-node=10 # 10 cores
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 
#SBATCH --mem 40G
#SBATCH -p standard
#SBATCH --account jcbnunez
#SBATCH --array=1-13

# This script is a pipeline which gather VCFs from all chromosomes.

#Load Modules
module load vcftools
module load bcftools
module load tabix

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/yey2sn/Overwintering_ms/16.Haplotypes
input_vcf=MC.DGRP.Taylor.merged.readyForImport.recode.vcf.gz

vcftools \
--gzvcf $input_vcf \
--from-bp 6671911 \
--to-bp 6671911 \
--chr 2L \
--recode --recode-INFO-all \
--out Spn27A

bgzip Spn27A.recode.vcf
tabix Spn27A.recode.vcf.gz

### DONE

echo "done" $(date)