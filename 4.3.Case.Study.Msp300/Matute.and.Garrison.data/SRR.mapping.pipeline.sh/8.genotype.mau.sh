#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=170G
#SBATCH --time=72:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez
#SBATCH --array=1-7

# This script will conduct genotype calling on the GenomeDBI object

#Load Modules
module load gatk

#Name of pipeline
PIPELINE=GenotypeGVCFs.mau

#Reference genome
REFERENCE=/project/berglandlab/Dmel_genomic_resources/Other_drosophilids/GCF_004382145.1_Dmauritiana_ASM438214v1_genomic.fna

#Intervals to analyze
intervals=mau.intervals.txt


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
###########################################################################
# Generate Folders and files
###########################################################################
###########################################################################

#Interval to analyze
i=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`

echo ${i} "is being processed" $(date)

# Identify the Genome database to genotyoe
GenomeDB_path=`echo /project/berglandlab/Drosophilids_datasets_yak_sech_sim_maur/Genomics_db_mau/DB_${i}`

echo "now processing DB" ${i} $(date)

###########################################################################
###########################################################################
# Genotype call the samples in the DBI merged set
###########################################################################
###########################################################################
mkdir ./out.vcfs

 gatk --java-options "-Xmx${JAVAMEM}" GenotypeGVCFs \
   -R $REFERENCE \
   -V gendb://$GenomeDB_path \
   -O ./out.vcfs/${i}.genotyped.raw.vcf.gz

echo ${i} "done" $(date)

