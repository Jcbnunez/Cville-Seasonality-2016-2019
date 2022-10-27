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
dgrp_vcf2l=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz

KAR=INVvsSTD
#File containing populations
std_samps=Std_samps_OnlyNames.txt
inv_samps=Inv_samps_OnlyNames.txt

dgrp_std=STD_DGRP_OnlyNames.txt
dgrp_inv=INV_DGRP_OnlyNames.txt
#samps=Het_samps_OnlyNames.txt

#window parameters
window=100000
step=50000

window2=500000
step2=100000

window3=10000
step3=5000

window4=5000
step4=1000

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

vcftools \
--gzvcf $input_vcf \
--weir-fst-pop $std_samps \
--weir-fst-pop $inv_samps \
--fst-window-size $window2 \
--fst-window-step $step2 \
--chr 2L \
--out FST.W_$window2.S_$step2.$KAR

vcftools \
--gzvcf $input_vcf \
--weir-fst-pop $std_samps \
--weir-fst-pop $inv_samps \
--fst-window-size $window3 \
--fst-window-step $step3 \
--chr 2L \
--out FST.W_$window3.S_$step3.$KAR

vcftools \
--gzvcf $input_vcf \
--weir-fst-pop $std_samps \
--weir-fst-pop $inv_samps \
--fst-window-size $window4 \
--fst-window-step $step4 \
--chr 2L \
--out FST.W_$window4.S_$step4.$KAR


#### same analysis in the DGRP
vcftools \
--gzvcf $dgrp_vcf2l \
--weir-fst-pop $dgrp_std \
--weir-fst-pop $dgrp_inv \
--fst-window-size $window \
--fst-window-step $step \
--chr 2L \
--out DGRP.FST.W_$window.S_$step.$KAR

#### same analysis in the DGRP -- wide window
vcftools \
--gzvcf $dgrp_vcf2l \
--weir-fst-pop $dgrp_std \
--weir-fst-pop $dgrp_inv \
--fst-window-size $window2 \
--fst-window-step $step2 \
--chr 2L \
--out DGRP.FST.W_$window2.S_$step2.$KAR

#### same analysis in the DGRP -- wide window
vcftools \
--gzvcf $dgrp_vcf2l \
--weir-fst-pop $dgrp_std \
--weir-fst-pop $dgrp_inv \
--fst-window-size $window3 \
--fst-window-step $step3 \
--chr 2L \
--out DGRP.FST.W_$window3.S_$step3.$KAR

vcftools \
--gzvcf $dgrp_vcf2l \
--weir-fst-pop $dgrp_std \
--weir-fst-pop $dgrp_inv \
--fst-window-size $window4 \
--fst-window-step $step4 \
--chr 2L \
--out DGRP.FST.W_$window4.S_$step4.$KAR


echo "done" $(date)