#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=8:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez

# This script will gather all sample data across the many lanes of sequencing

#Load Rivanna modules 
module load samtools
module load qualimap

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/yey2sn/Alyssa_single_inds/fastq_to_VCF

#Define important file locations
JOINT_BAMS=/scratch/yey2sn/Alyssa_single_inds/fastq_to_VCF/joint_bams

#Name of pipeline
PIPELINE=Merge_Bams

# This is a file with the name all the samples to be processed. one sample name per line with all the infor
#this is a file has unique file identifiers
SAMPLE_FILE=/scratch/yey2sn/Alyssa_single_inds/fastq_to_VCF/Alyssa_ind_logic_comparison_inds.txt

#						File.1 Files                      File.2  Files Logical		Unique names
#CM.002.0916_SL256659.fastq.gz	8	CM.002.0916_SL256659.fastq.gz	8	yes	yes8	CM.002.0916_SL256659
#CM.003.0722_SL256747.fastq.gz	8	CM.003.0722_SL256747.fastq.gz	8	yes	yes8	CM.003.0722_SL256747
#CM.004.0819_SL256728.fastq.gz	8	CM.004.0819_SL256728.fastq.gz	8	yes	yes8	CM.004.0819_SL256728

####################
#Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
QUAL=40 # Quality threshold for samtools
JAVAMEM=18g # Java memory
####################

###########################################################################
###########################################################################
# Determine sample to process, "i" and read files
###########################################################################
###########################################################################
 
 #Column 7 contains the unique name
 
i=`awk -F "\t" '{print $7}' $SAMPLE_FILE | sed -n ${SLURM_ARRAY_TASK_ID}p`

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

if [ -d "Merged_Bams" ]
then
	echo "Working mapping_stats folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/Merged_Bams
	date
fi

if [ -d "Merged_Bams_qualimap" ]
then
	echo "Working mapping_stats folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/Merged_Bams_qualimap
	date
fi


###########################################################################
# Start pipeline
###########################################################################
# Lets do some light trimming of the reads
###########################################################################

echo "I will merge these files"  $JOINT_BAMS/*_${i}.joint.srt.rmdp.bam

#Make temporary linefile
ls $JOINT_BAMS/*_${i}.joint.srt.rmdp.bam > ${i}.guide.txt

samtools merge \
-b ${i}.guide.txt \
$WORKING_FOLDER/Merged_Bams/${i}.Lanes_merged.bam

#remove the temporary guide file
rm ${i}.guide.txt

# Assess quality of final file
qualimap bamqc \
 -bam $WORKING_FOLDER/Merged_Bams/${i}.Lanes_merged.bam   \
 -outdir $WORKING_FOLDER/Merged_Bams_qualimap/Qualimap_LaneMerged_${i} \
 --java-mem-size=$JAVAMEM

###########################################################################
###########################################################################
# Inform that sample is done
###########################################################################
###########################################################################
# This part of the pipeline will produce a notification of completion. 

echo ${i} " completed" >> $WORKING_FOLDER/${PIPELINE}.completion.log

echo "pipeline completed" $(date)
