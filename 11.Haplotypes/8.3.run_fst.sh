#! /bin/bash

#SBATCH -J fst_inv_std 
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
#dgrp_vcf2l=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz

KAR=INVvsSTD
#File containing populations
std_samps=Std_samps_OnlyNames.txt
inv_samps=Inv_samps_OnlyNames.txt
#samps=Het_samps_OnlyNames.txt

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

vcftools \
--gzvcf $input_vcf \
--weir-fst-pop $std_samps \
--weir-fst-pop $inv_samps \
--fst-window-size $window \
--fst-window-step $step \
--chr 2L \
--out FST.W_$window.S_$step.$KAR

sed -i "s/^/${i}\t/" FST.W_$window.S_$step.$KAR.windowed.pi

echo "done" $(date)