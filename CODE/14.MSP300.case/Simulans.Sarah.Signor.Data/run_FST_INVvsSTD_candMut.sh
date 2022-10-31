#! /bin/bash

#SBATCH -J fst.sim 
#SBATCH --ntasks-per-node=10 # 10 cores
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 
#SBATCH --mem 40G
#SBATCH -p standard
#SBATCH --account jcbnunez


# This script is a pipeline which gather VCFs from all chromosomes.

#Load Modules
module load vcftools


##guide=8.vcftools_guide_dat.txt
###

input_vcf=Dsim_all_2L.maf0.01.recode.vcf.gz

###########################################################################
###########################################################################
# Begin Pipeline
###########################################################################
###########################################################################
#This part of the pipeline will generate log files to record warnings and completion status

###########################################################################
###########################################################################
# Filter final VCF
###########################################################################
###########################################################################

vcftools --gzvcf $input_vcf \
--weir-fst-pop homozygous.2.txt \
--weir-fst-pop homozygous.0.txt \
--out ./sim.mut.haps.fst \
--fst-window-size 10000 \
--fst-window-step 5000 \
--maf 0.01 

#####
### DONE

echo "done" $(date)