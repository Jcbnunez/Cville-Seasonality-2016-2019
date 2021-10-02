#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=170G
#SBATCH --time=72:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez

# This script will conduct genotype calling on the GenomeDBI object

#Load Modules
module load gatk

#Name of pipeline
PIPELINE=GenotypeGVCFs

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/yey2sn/phasing_droso_data

#Reference genome
REFERENCE=/project/berglandlab/Dmel_fasta_refences/holo_dmel_6.12.fa

#Intervals to analyze
intervals=/scratch/yey2sn/phasing_droso_data/Intervals_Dmel.txt


###########################################################################
#Parameters
#Java
JAVAMEM=160G
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
# Generate Folders and files
###########################################################################
###########################################################################

#Interval to analyze
i=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`

echo ${i} "is being processed" $(date)

# Identify the Genome database to genotyoe
GenomeDB_path=`echo /project/berglandlab/Dmel_Single_Individuals/Genotype_Database_GATK/DB_${i}`

echo "now processing DB" ${i} $(date)

###########################################################################
###########################################################################
# Genotype call the samples in the DBI merged set
###########################################################################
###########################################################################

 gatk --java-options "-Xmx${JAVAMEM}" GenotypeGVCFs \
   -R $REFERENCE \
   -V gendb://$GenomeDB_path \
   -O $WORKING_FOLDER/${i}.genotyped.raw.vcf.gz

echo ${i} "done" $(date)