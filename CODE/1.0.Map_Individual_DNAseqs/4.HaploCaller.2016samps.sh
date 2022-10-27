#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez

# This script will initiate a pipeline which will add read group info and index bams. It will then proceed to call haplotypes (gVCFs)
# Prepared by Joaquin C. B. Nunez, PhD -- Sep 25, 2020 + updated April 1, 2021
# yey2sn@virginia.edu

#Load Modules
module load gatk
module load picard
module load tabix

#Name of pipeline
PIPELINE=Haplocaller

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/yey2sn/Alyssa_single_inds/fastq_to_VCF

# User defined inputs -- this represents the name of the samples
BIO_SAMPLES=/scratch/yey2sn/Alyssa_single_inds/fastq_to_VCF/Alyssa_ind_logic_comparison_inds.txt

#Where the bam files are located
BAMS_FOLDER=/scratch/yey2sn/Alyssa_single_inds/fastq_to_VCF/Merged_Bams

#Reference genome
REFERENCE=/project/berglandlab/Dmel_fasta_refences/holo_dmel_6.12.fa

#Sample suffixes and post-fixes. What tags are expected across all samples?
# Understanding of this comes from the previous pipeline
SUFFIX="Lanes_merged"

#Parameters

#Java
JAVAMEM=18G
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE

#Read Information
Group_library="Alyssa_2016"
Library_Platform="illumina"
Group_platform="A2016"

#HaploCaller -- heterozygocity prior
HET=0.005

###########################################################################
###########################################################################
# Determine sample to process, "i"
###########################################################################
###########################################################################
 
i=`awk -F "\t" '{print $7}' $BIO_SAMPLES | sed -n ${SLURM_ARRAY_TASK_ID}p`

###########################################################################
###########################################################################
# Begin Pipeline
###########################################################################
###########################################################################
#This part of the pipeline will generate log files to record warnings and completion status


# Welcome message
echo "your unique run id is" $unique_run_id

if [[ -e "${PIPELINE}.completion.log" ]]
then
	echo "Warning log exist"
	echo "lets move on"
	date
else 
	echo "Log doesnt exist. lets fix that"
	touch $WORKING_FOLDER/${PIPELINE}.completion.log
	date
fi

# Move to working directory
cd $WORKING_FOLDER

###########################################################################
###########################################################################
# Generate Folders and files
###########################################################################
###########################################################################
# this part of the script will check and generate, if necessary, all of the output folders used in the script

echo "have you checked if the folders where already built with mkdir?"

if [[ -d "RGSM_final_bams" ]]
then
	echo "Working RGSM_final_bams folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/RGSM_final_bams
	date
fi


if [[ -d "haplotype_calling" ]]
then
	echo "Working haplotype_calling folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/haplotype_calling
	date
fi

###########################################################################
###########################################################################
# Forcing a uniform read group to the joint bam file
###########################################################################
###########################################################################

java -jar $PICARD AddOrReplaceReadGroups \
	I=$BAMS_FOLDER/${i}.$SUFFIX.bam \
	O=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam \
	RGLB=$Group_library \
	RGPL=$Library_Platform \
	RGPU=$Group_platform \
	RGSM=${i}

###########################################################################
###########################################################################
# Index Bam files
###########################################################################
###########################################################################

java -jar $PICARD BuildBamIndex \
      I=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam

###########################################################################
###########################################################################
# Haplotype Calling
###########################################################################
###########################################################################
# Call haplotypes with GATK

gatk --java-options "-Xmx${JAVAMEM}" HaplotypeCaller \
	-R $REFERENCE \
	-I $WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam \
	-O $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf \
	--heterozygosity $HET \
	-ERC GVCF 

###########################################################################
###########################################################################
# Compress and index with Tabix
###########################################################################
###########################################################################

bgzip $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf
tabix $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf.gz

echo ${i} "completed" $(date) >> $WORKING_FOLDER/${PIPELINE}.completion.log

echo "done" $(date)
