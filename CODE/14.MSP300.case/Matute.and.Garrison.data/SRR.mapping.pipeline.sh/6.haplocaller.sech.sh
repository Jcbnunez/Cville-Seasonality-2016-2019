#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez
#SBATCH -o ./slurmOut/map%A_%a.out # Standard output
#SBATCH -e ./slurmOut/map%A_%a.err # Standard error
#SBATCH --array=1-41


# This script will initiate a pipeline which will add read group info and index bams. It will then proceed to call haplotypes (gVCFs)
# Prepared by Joaquin C. B. Nunez, PhD -- Sep 25, 2020 + updated April 1, 2021
# yey2sn@virginia.edu

#Load Modules
module load gatk
module load picard
module load tabix

#Sample suffixes and post-fixes. What tags are expected across all samples?
# Understanding of this comes from the previous pipeline

#Parameters

#Java
JAVAMEM=18G
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE

#HaploCaller -- heterozygocity prior
HET=0.005

###########################################################################
###########################################################################
# Determine sample to process, "i"
###########################################################################
###########################################################################
 
REFERENCE=/project/berglandlab/Dmel_genomic_resources/Other_drosophilids/GCF_004382195.2_Dsechelia_ASM438219v2_genomic.fna
###
Metadat=sech.metadat.mapping.txt
OUT=/scratch/yey2sn/Overwintering_ms/msp300.case/mapping_output


###
sampleId=`sed '1d' $Metadat | awk -F "\t" '{print $1}' | sed -n ${SLURM_ARRAY_TASK_ID}p`
config=`sed '1d' $Metadat | awk -F "\t" '{print $2}' | sed -n ${SLURM_ARRAY_TASK_ID}p`
sp=`sed '1d' $Metadat | awk -F "\t" '{print $3}' | sed -n ${SLURM_ARRAY_TASK_ID}p`
SRR=`sed '1d' $Metadat | awk -F "\t" '{print $4}' | sed -n ${SLURM_ARRAY_TASK_ID}p`
ReadF=`sed '1d' $Metadat | awk -F "\t" '{print $5}' | sed -n ${SLURM_ARRAY_TASK_ID}p`
ReadR=`sed '1d' $Metadat | awk -F "\t" '{print $6}' | sed -n ${SLURM_ARRAY_TASK_ID}p`

#Read Information
Group_library="sech_lib"
Library_Platform="illumina"
Group_platform="sech_group"


###########################################################################
###########################################################################
# Forcing a uniform read group to the joint bam file
###########################################################################
###########################################################################

java -jar $PICARD AddOrReplaceReadGroups \
	I=$OUT/$sp/$SRR/$sampleId.srt.rmdp.bam \
	O=$OUT/$sp/$SRR/$sampleId.srt.rmdp.RG.bam \
	RGLB=$Group_library \
	RGPL=$Library_Platform \
	RGPU=$Group_platform \
	RGSM=$sampleId

###########################################################################
###########################################################################
# Index Bam files
###########################################################################
###########################################################################

java -jar $PICARD BuildBamIndex \
      I=$OUT/$sp/$SRR/$sampleId.srt.rmdp.RG.bam

###########################################################################
###########################################################################
# Haplotype Calling
###########################################################################
###########################################################################
# Call haplotypes with GATK

gatk --java-options "-Xmx${JAVAMEM}" HaplotypeCaller \
	-R $REFERENCE \
	-I $OUT/$sp/$SRR/$sampleId.srt.rmdp.RG.bam \
	-O $OUT/$sp/$SRR/$sampleId.raw.g.vcf \
	--heterozygosity $HET \
	-ERC GVCF 

###########################################################################
###########################################################################
# Compress and index with Tabix
###########################################################################
###########################################################################

bgzip $OUT/$sp/$SRR/$sampleId.raw.g.vcf
tabix $OUT/$sp/$SRR/$sampleId.raw.g.vcf.gz

echo "done" $(date)
