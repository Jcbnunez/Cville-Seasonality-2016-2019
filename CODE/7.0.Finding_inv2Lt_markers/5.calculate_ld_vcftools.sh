#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=40G
#SBATCH --time=8:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab

###########################################################################
###########################################################################
# Load R RIVANNA modules
###########################################################################
###########################################################################

#Load Modules
module load vcftools
module load bcftools
module load tabix

###########################################################################
###########################################################################
# Internal variables
###########################################################################
###########################################################################

WORKING_FOLDER=/project/berglandlab/Dmel_genomic_resources/Inversions/DGRP_2lt_Markers
#Final VCF file
IN_GZVCF=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz
#File containing populations
samples_to_keep=./lines_used_for_2lt_demarcation.txt
#determine id
id=Rank99rho
#SNP ids
position_file=./SNPids_99rank.txt

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
# Select Population file
###########################################################################
###########################################################################

#i=`sed -n ${SLURM_ARRAY_TASK_ID}p $population_guide`

###########################################################################
###########################################################################
# Filter final VCF
###########################################################################
###########################################################################

vcftools \
--gzvcf $IN_GZVCF \
--keep $samples_to_keep \
--geno-r2-positions $position_file \
--out $id

echo "done" $(date)

