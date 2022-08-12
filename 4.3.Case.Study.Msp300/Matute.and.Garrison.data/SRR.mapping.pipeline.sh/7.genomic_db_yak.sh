#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=176G
#SBATCH --time=72:00:00
#SBATCH --partition=largemem
#SBATCH --account=jcbnunez
#SBATCH --array=1-7

# This script will merge gVCFs into a unified database for genotype calling.
# Because this and following steps are so memory intensive, this will be done using
# a per chromosome approach

#Load Modules
module load gatk

#Name of pipeline
PIPELINE=GenomicsDBImport.yak

# User defined inputs -- this represents the name of the samples
sp=yak
sample_map=/scratch/yey2sn/Overwintering_ms/msp300.case/yak_samps.to.haplo.txt
#This file looks like this
#  sample1      sample1.vcf.gz
#  sample2      sample2.vcf.gz
#  sample3      sample3.vcf.gzs
#Intervals to analyze for each species
intervals=yak.intervals.txt

######
#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/yey2sn/Overwintering_ms/msp300.case/
DB_location=/project/berglandlab/Drosophilids_datasets_yak_sech_sim_maur/Genomics_db_yak

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
       --tmp-dir $WORKING_FOLDER/TEMP_MERGEVCF_${i} \
       --reader-threads $CPU \
       -L ${i}

echo ${i} "done" $(date)
