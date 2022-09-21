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

guide=8.vcftools_guide_dat.txt
###

input_vcf=$( cat $guide  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $1 }' )
samps=$( cat $guide  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $2 }' )
window=$( cat $guide  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $3 }' )
step=$( cat $guide  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $4 }' )
KAR=$( cat $guide  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $5 }' )
pop=$( cat $guide  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $6 }' )

### Bed file for hap count
BED=bed_guide.file.2L.BED

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
--maf 0.05 \
--keep $samps \
--hapcount $BED \
--chr 2L \
--out HAPCOUNT.$pop.W_$window.S_$step.$KAR

### DONE

echo "done" $(date)