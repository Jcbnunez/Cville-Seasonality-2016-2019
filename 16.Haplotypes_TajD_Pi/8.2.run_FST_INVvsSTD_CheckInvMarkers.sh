#! /bin/bash

#SBATCH -J Std_pi_d 
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

##guide=8.vcftools_guide_dat.txt
###

input_vcf=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz

PosInt=./SNPs.for.fst.analysis.ids.only.txt
head $PosInt
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

vcftools --gzvcf $input_vcf \
--weir-fst-pop CM.STD.samps.txt \
--weir-fst-pop CM.INV.samps.txt \
--out ./MARKERS.inv.vs.std.cm.fst \
--maf 0.01 \
--positions $PosInt



### DONE

echo "done" $(date)