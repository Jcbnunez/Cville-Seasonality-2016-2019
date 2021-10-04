#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=8:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez


# This script will initiate a pipeline which will do some quality QC on the reads and then will proceed to map the reads to a reference genome.
# Prepared by Joaquin C. B. Nunez, PhD -- March 29, 2021
# yey2sn@virginia.edu

# NOTES ON NOMENCLATURE: This script uses nomenclature which can be confusing. The first part of the script split raw reads into insert-"merged"-reads (hereby called merged) and unmerged reads (those which were not merged). As such, all operations done using ether of these reads will have the term "merged" or "unmerged" attached to them. At a later point in the script, I combine bam files using "samtools merge" the output of this combination is a joint-bam file (hereby called "joint"). Thus, the joint suffix referes to this step. Other suffix used here are: "srt" which mean picard sorted, and "rmdp" which mean picard-removed duplicated reads.
  
#Load Rivanna modules 
module load gcc/9.2.0
module load bwa/0.7.17
module load bbmap
module load fastqc
module load samtools
module load qualimap
module load picard

#Define important file locations
#RAW READS indicates the folder where the raw reads are stored.
RAW_READS=/project/berglandlab/alyssa/wildDmel2016/fastq

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/yey2sn/Alyssa_single_inds/fastq_to_VCF

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/project/berglandlab/Dmel_fasta_refences/holo_dmel_6.12.fa

#This is a unique number id which identifies this run
unique_run_id=`date +%N`

#Name of pipeline
PIPELINE=Merge_and_Trim

####################
#Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
QUAL=40 # Quality threshold for samtools
JAVAMEM=18g # Java memory
####################

## PREPARE GUIDE FILES
## Read guide files
# This is a file with the name all the samples to be processed. one sample name per line with all the infor
SAMPLE_FILE=/scratch/yey2sn/Alyssa_single_inds/fastq_to_VCF/Alyssa_ind_reads_guideFile.txt

#Example: -- the headers are just examples. The actual file has no header
 ##                                          File1       lane_1                      sample_1
## 5  HL3CVALXX_s7_1_CM.002.0916_SL256659.fastq.gz HL3CVALXX_s7 CM.002.0916_SL256659.fastq.gz
## 7  HL3CVALXX_s7_1_CM.003.0722_SL256747.fastq.gz HL3CVALXX_s7 CM.003.0722_SL256747.fastq.gz
## 10 HL3CVALXX_s7_1_CM.004.0819_SL256728.fastq.gz HL3CVALXX_s7 CM.004.0819_SL256728.fastq.gz

##                                           File2       lane_2                      sample_2
## 5  HL3CVALXX_s7_2_CM.002.0916_SL256659.fastq.gz HL3CVALXX_s7 CM.002.0916_SL256659.fastq.gz
## 7  HL3CVALXX_s7_2_CM.003.0722_SL256747.fastq.gz HL3CVALXX_s7 CM.003.0722_SL256747.fastq.gz
## 10 HL3CVALXX_s7_2_CM.004.0819_SL256728.fastq.gz HL3CVALXX_s7 CM.004.0819_SL256728.fastq.gz

##    sample_name      SL_extension                       Merged_name
## 5  CM.002.0916 SL256659.fastq.gz HL3CVALXX_s7_CM.002.0916_SL256659
## 7  CM.003.0722 SL256747.fastq.gz HL3CVALXX_s7_CM.003.0722_SL256747
## 10 CM.004.0819 SL256728.fastq.gz HL3CVALXX_s7_CM.004.0819_SL256728

###########################################################################
###########################################################################
# Determine sample to process, "i" and read files
###########################################################################
###########################################################################
 
i=`awk -F "\t" '{print $9}' $SAMPLE_FILE | sed -n ${SLURM_ARRAY_TASK_ID}p`
read1=`awk -F "\t" '{print $1}' $SAMPLE_FILE | sed -n ${SLURM_ARRAY_TASK_ID}p`
read2=`awk -F "\t" '{print $4}' $SAMPLE_FILE | sed -n ${SLURM_ARRAY_TASK_ID}p`

###########################################################################
###########################################################################
# Begin Pipeline
###########################################################################
###########################################################################
#This part of the pipeline will generate log files to record warnings and completion status

# Welcome message
echo "your unique run id is" $unique_run_id

if [[ -e "${PIPELINE}.warnings.log" ]]
then
	echo "Warning log exist"
	echo "lets move on"
	date
else 
	echo "Log doesnt exist. lets fix that"
	touch $WORKING_FOLDER/${PIPELINE}.warnings.log
	date
fi

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

#Generating new folders 
echo "have you checked if the folders where already built with mkdir?"
if [[ -d "merged_reads" ]]
then
	echo "Working merged_reads folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/merged_reads
	date
fi

if [ -d "unmerged_reads" ]
then
	echo "Working unmerged_reads folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/unmerged_reads
	date
fi

if [ -d "fastqc_merged" ]
then
	echo "Working unmerged_reads folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/fastqc_merged
	date
fi

###########################################################################
###########################################################################
# Trim and merge reads
###########################################################################
###########################################################################
# This part of the pipeline will trim and merge the reads. It is very likely that the reads will be split into merged and unmerged. Both reads will be mapped. This loop operates using a while-read-do-done structure. the while loop is feed a file "SAMPLE_FILE" where  all sample names are stored, one name per line. This can be leveraged for parallelization.

	echo ${i} "is now processing"
	date

	mkdir $WORKING_FOLDER/merged_reads/${i}
	mkdir $WORKING_FOLDER/unmerged_reads/${i}
	
	echo "now merging reads for" ${i}
	
	bbmerge.sh \
	in1=${RAW_READS}/$read1 in2=${RAW_READS}/$read2 \
	out=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq \
	outu1=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.1.fq \
	outu2=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.2.fq \
	-strict
	
		#Sanity checks	
	if [ -s $WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq ]; then
	echo ${i} "merged reads file is not empty... thats good"
	else
	echo "File is empty -- WARNING ISSUED!"
	echo ${i} "Merged reads is empty! check the script, debug, and rerun" >> $WORKING_FOLDER/${PIPELINE}.warnings.log
	fi
	
	if [ -s $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.1.fq ]; then
	echo ${i} "Pair 1 reads file is not empty... thats good"
	else
	echo "File is empty -- WARNING ISSUED!"
	echo ${i} "Pair 1 reads is empty! check the script, debug, and rerun" >> $WORKING_FOLDER/${PIPELINE}.warnings.log
	fi
	
	if [ -s $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.2.fq ]; then
	echo ${i} "Pair 2 reads file is not empty... thats good"
	else
	echo "File is empty -- WARNING ISSUED!"
	echo ${i} "Pair 2 reads is empty! check the script, debug, and rerun" >> $WORKING_FOLDER/${PIPELINE}.warnings.log
	fi	

############## Now lest do some QC on the reads

	fastqc $WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq \
	--outdir $WORKING_FOLDER/fastqc_merged

###########################################################################
###########################################################################
# Inform that sample is done
###########################################################################
###########################################################################
# This part of the pipeline will produce a notification of completion. 

echo ${i} " completed" >> $WORKING_FOLDER/${PIPELINE}.completion.log

echo "pipeline" ${PIPELINE} $(date)
