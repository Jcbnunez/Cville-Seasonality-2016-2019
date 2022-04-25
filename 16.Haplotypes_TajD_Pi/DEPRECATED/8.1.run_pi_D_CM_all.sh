#! /bin/bash

#SBATCH -J vcftools_all 
#SBATCH --ntasks-per-node=10 # 10 cores
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 
#SBATCH --mem 40G
#SBATCH -p standard
#SBATCH --account jcbnunez



# This script is a pipeline which gather VCFs from all chromosomes.

#Load Modules
module load vcftools
module load bcftools
module load tabix

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/yey2sn/Overwintering_ms/11.Haplotypes

#Final VCF file
input_vcf=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz
dgrp_vcf2l=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz

KAR=All

#window parameters
window=100000
step=50000

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

###########################################################################
###########################################################################
# Filter final VCF
###########################################################################
###########################################################################

vcftools \
--gzvcf $input_vcf \
--window-pi $window \
--window-pi-step $step \
--chr 2L \
--out Pi.CM.W_$window.S_$step.$KAR


vcftools \
--gzvcf $input_vcf \
--TajimaD $window \
--chr 2L \
--out D.CM.W_$window.S_$step.$KAR

##################
#### DGRP SAMPLEs
##################

vcftools \
--gzvcf $dgrp_vcf2l \
--window-pi $window \
--window-pi-step $step \
--chr 2L \
--out Pi.DGRP.W_$window.S_$step.$KAR

vcftools \
--gzvcf $dgrp_vcf2l \
--TajimaD $window \
--chr 2L \
--out D.DGRP.W_$window.S_$step.$KAR

### DONE

echo "done" $(date)