#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=8:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez

# This script is the second step in the qc-trim-map of the Alyssa's data

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
PIPELINE=Trim_and_Map

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

###########################################################################
###########################################################################


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

if [ -d "mapping_stats" ]
then
	echo "Working mapping_stats folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/mapping_stats
	date
fi

if [ -d "read_stats" ]
then
	echo "Working read_stats folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/read_stats
	date
fi

if [ -d "joint_bams" ]
then
	echo "Working joint_bams folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/joint_bams
	date
fi

if [ -d "joint_bams_qualimap" ]
then
	echo "Working joint_bams_qualimap folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/joint_bams_qualimap
	date
fi

###########################################################################
# Start pipeline
###########################################################################
# Lets do some light trimming of the reads
###########################################################################

### Work on merged reads
### Decide on the trimming parameters based on fastQC step done before this script.

	echo ${i} "Trimming merged reads"

	bbduk.sh \
	in=`echo $WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq` \
	out=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.trim.fq \
	ftl=15 ftr=285 qtrim=w trimq=20

	rm  $WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq

### Work on unmerged reads
	echo ${i} "Trimming unmerged reads"
	
	bbduk.sh \
	in=`echo $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.1.fq` \
	in2=`echo $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.2.fq` \
	out=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.trim.1.fq \
	out2=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.trim.2.fq \
	ftl=15 qtrim=w trimq=20
	
	rm $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.1.fq
	rm $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.2.fq

###########################################################################
###########################################################################
# Map reads to a reference
###########################################################################
###########################################################################
# this part will map reads to the reference genome. Because the reads are likely split into two groups, this script will loop over both types of reads. After reads have been mapped, they will be compressed into bam files, sorted, and duplicates will be removed. I will also conduct an intermediary QC step with Qualimap. Because there are inherent QC steps here, I have avoided adding extra "warnings" in the log. Remember to take a look at the qualimap and the flagstat outputs to check for inconsistencies.  

for j in merged unmerged
	do # Begin loop of j
	
		########################################
#J loop#	# Starting mapping
	echo "I will first map ${j} reads of" ${i}
	
#J loop#	# I will conduct the mapping with BWA-MEM
	
	if [[ ${j} == "merged" ]]; then
		echo "seems this is merged data, lets map it"
		bwa mem \
		-M \
		-t $CPU \
		$REFERENCE \
		$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.reads.strict.trim.fq \
		> $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.sam
		
	elif [[ ${j} == "unmerged" ]]; then
		echo "seems this is unmerged data, lets map it using a 1-2 approach"
		bwa mem \
		-M \
		-t $CPU \
		$REFERENCE \
		$WORKING_FOLDER/unmerged_reads/${i}/${i}.${j}.reads.trim.1.fq \
		$WORKING_FOLDER/unmerged_reads/${i}/${i}.${j}.reads.trim.2.fq \
		> $WORKING_FOLDER/unmerged_reads/${i}/${i}.${j}.sam
	
	else
		echo "I cant tell what type of data is this -- WARNING!"
		echo ${i} "Something is wrong at the mapping stage" $(date) \
		  $Project_name.warnings.$unique_run_id.log
	fi

#J loop#	#I will now extract some summary stats
	samtools flagstat \
	--threads $CPU \
	$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.sam \
	> $WORKING_FOLDER/${j}_reads/${i}/${i}.flagstats_raw_${j}.sam.txt

#J loop#	#build bam files
	samtools view \
	-b \
	-q $QUAL \
	--threads $CPU  \
	$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.sam \
	> $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.bam

#J loop#	# Sort with picard
	# Notice that once a file has been sorted it is added the "srt" suffix
	java -Xmx$JAVAMEM \
	-jar $PICARD SortSam \
	I=$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.bam \
	O=$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.bam \
	SO=coordinate \
	VALIDATION_STRINGENCY=SILENT

#J loop# Remove duplicates with picard
	# Notice that once a file has been sorted it is added the "rmdp" suffix
	java -Xmx$JAVAMEM \
	-jar $PICARD MarkDuplicates \
	I=$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.bam \
	O=$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.rmdp.bam \
	M=$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.dupstat.txt \
	VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

#J loop# Lets do QC on the bam file
	qualimap bamqc \
	-bam $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.rmdp.bam  \
	-outdir $WORKING_FOLDER/mapping_stats/Qualimap_${i} \
	--java-mem-size=$JAVAMEM

#J loop#	# Clean intermediate files
	rm $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.sam
	rm $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.bam
	rm $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.bam

#J loop#	# Housekeeping
	mv $WORKING_FOLDER/${j}_reads/${i}/${i}.flagstats_raw_${j}.sam.txt \
	$WORKING_FOLDER/mapping_stats
	mv $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.dupstat.txt \
	$WORKING_FOLDER/mapping_stats

	mv $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.reads.strict.trim_fastqc.html \
	$WORKING_FOLDER/read_stats
	mv $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.reads.strict.trim_fastqc.zip \
	$WORKING_FOLDER/read_stats

#J loop#	
	done # End loop of j

###########################################################################
###########################################################################
# Merge and asses the final file
###########################################################################
###########################################################################
# Here I will merge the bam outputs from the merge and unmerged portions of the pipeline. Subsequently, I will once again sort and remove duplicated, before performing the final QC on the aligment.

# Merge bams
java -Xmx$JAVAMEM \
 -jar $PICARD MergeSamFiles \
 I=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.srt.rmdp.bam  \
 I=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.srt.rmdp.bam  \
 O=$WORKING_FOLDER/joint_bams/${i}.joint.bam

# Sort merge bams
java -Xmx$JAVAMEM \
 -jar $PICARD SortSam \
 I=$WORKING_FOLDER/joint_bams/${i}.joint.bam \
 O=$WORKING_FOLDER/joint_bams/${i}.joint.srt.bam \
 SO=coordinate \
 VALIDATION_STRINGENCY=SILENT

# Remove duplicates of final file
java -Xmx$JAVAMEM \
 -jar $PICARD MarkDuplicates \
 I=$WORKING_FOLDER/joint_bams/${i}.joint.srt.bam \
 O=$WORKING_FOLDER/joint_bams/${i}.joint.srt.rmdp.bam  \
 M=$WORKING_FOLDER/mapping_stats/${i}.joint.dupstat.txt \
 VALIDATION_STRINGENCY=SILENT \
 REMOVE_DUPLICATES=true

# Assess quality of final file
qualimap bamqc \
 -bam $WORKING_FOLDER/joint_bams/${i}.joint.srt.rmdp.bam   \
 -outdir $WORKING_FOLDER/joint_bams_qualimap/Qualimap_JointBam_${i} \
 --java-mem-size=$JAVAMEM
 
# Remove intermediary files
rm $WORKING_FOLDER/joint_bams/${i}.joint.bam
rm $WORKING_FOLDER/joint_bams/${i}.joint.srt.bam

###########################################################################
###########################################################################
# Inform that sample is done
###########################################################################
###########################################################################
# This part of the pipeline will produce a notification of completion. 

echo ${i} " completed" >> $WORKING_FOLDER/${PIPELINE}.completion.log

echo "pipeline completed" $(date)
