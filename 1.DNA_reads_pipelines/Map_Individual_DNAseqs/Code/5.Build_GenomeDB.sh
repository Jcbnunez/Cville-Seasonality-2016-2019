#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=170G
#SBATCH --time=72:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez

# This script will merge gVCFs into a unified database for genotype calling.
# Because this and following steps are so memory intensive, this will be done using
# a per chromosome approach

#Load Modules
module load gatk

#Name of pipeline
PIPELINE=GenomicsDBImport

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/yey2sn/phasing_droso_data

# User defined inputs -- this represents the name of the samples
sample_map=/scratch/yey2sn/phasing_droso_data/samples_to_haplotype.txt

#This file looks like this
#  sample1      sample1.vcf.gz
#  sample2      sample2.vcf.gz
#  sample3      sample3.vcf.gzs

#Intervals to analyze
intervals=/scratch/yey2sn/phasing_droso_data/Intervals_Dmel.txt

DB_location=/project/berglandlab/Dmel_Single_Individuals/Genotype_Database_GATK

###########################################################################
#Parameters
#Java
JAVAMEM=150G
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
###########################################################################


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
# check samples
###########################################################################
###########################################################################

for k in $(seq 1 $( wc -l < $sample_map ))
do

name=`awk -F "\t" '{print $1}' $sample_map | sed -n ${k}p`
file=`awk -F "\t" '{print $2}' $sample_map | sed -n ${k}p`


if [[ -e "${file}" ]]
then
	status="exist"
else 
	status="WARNING"
fi

echo $name $status

done

###########################################################################
###########################################################################
# Generate Folders and files
###########################################################################
###########################################################################

#Interval to analyze
i=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`

echo ${i} "is being processed" $(date)

###########################################################################
###########################################################################
# Merge VCFs using GenomicsDBImport
###########################################################################
###########################################################################

mkdir $WORKING_FOLDER/TEMP_MERGEVCF_${i}

  gatk --java-options "-Xmx${JAVAMEM} -Xms${JAVAMEM}" \
       GenomicsDBImport \
       --genomicsdb-workspace-path $DB_location/DB_${i} \
       --batch-size 50 \
       --sample-name-map $sample_map \
       --tmp-dir=$WORKING_FOLDER/TEMP_MERGEVCF_${i} \
       --reader-threads $CPU \
       -L ${i}

echo ${i} "done" $(date)
