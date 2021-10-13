#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab

# This script is a pipeline which gather VCFs from all chromosomes.

#Load Modules
module load vcftools
module load bcftools
module load tabix

#Name of pipeline
PIPELINE=Inv2Lt_dat

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/project/berglandlab/Dmel_genomic_resources/Inversions/DGRP_2lt_Markers

#in-file VCFS
IN_GZVCF=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz

###########################################################################
###########################################################################
# Begin Pipeline
###########################################################################
###########################################################################
#This part of the pipeline will generate log files to record warnings and completion status

# Move to working directory
cd $WORKING_FOLDER

###########################################################################
###########################################################################
# Filter final VCF
###########################################################################
###########################################################################

FILTER=chr2L.miss80

vcftools \
--gzvcf $IN_GZVCF \
--chr 2L \
--remove-filtered-all \
--max-missing 0.8 \
--max-alleles 2 \
--recode \
--recode-INFO-all \
--out $PIPELINE.$FILTER

#bgzip and tabix
bgzip $PIPELINE.$FILTER.recode.vcf
tabix $PIPELINE.$FILTER.recode.vcf.gz
	
bcftools query $PIPELINE.$FILTER.recode.vcf.gz -l snp.vcf > $PIPELINE.$FILTER.sample_names.txt

echo "done" $(date)
