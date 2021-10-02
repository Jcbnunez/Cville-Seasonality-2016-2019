#!/usr/bin/env bash
#
#
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=85G
#SBATCH --time=36:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez

###########################################################################
#Parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
###########################################################################

### Note: where is Whatshap?
#Pip installs WhatsHap into $HOME/.local/bin. Then add $HOME/.local/bin to your $PATH and run the tool:

#Adding whapshapp to our path
export PATH=$HOME/.local/bin:$PATH

## Other rivanna dependencies
module load picard
module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822
export PATH=$HOME/.local/bin:$PATH

reference=/project/berglandlab/Dmel_fasta_refences/holo_dmel_6.12.fa

guide_file=/scratch/yey2sn/phasing_droso_data/Whatshap/samples_to_phase_chr.txt

PIPELINE=whapshapp_phase

### Master VCF file
MASTERVCF=/project/berglandlab/Dmel_Single_Individuals/Unphased_VCF_files/Autosomes_and_X_calibrated/Dmel_inds_AB2016_PE2018.vcf.gz

### Where the output will be saved
WORK_FOLDER=/scratch/yey2sn/phasing_droso_data/Whatshap/tmp

PHASED_FOLDER=/scratch/yey2sn/phasing_droso_data/Whatshap/pahsed_output

###############################

### generate completion log
if [[ -e "${PIPELINE}.completion.log" ]]
then
	echo "Warning log exist"
	echo "lets move on"
	date
else 
	echo "Log doesnt exist. lets fix that"
	touch ./${PIPELINE}.completion.log
	date
fi

### get parameters
sample=`awk -F "\t" '{print $1}' $guide_file | sed -n ${SLURM_ARRAY_TASK_ID}p`
bamfile=`awk -F "\t" '{print $2}' $guide_file | sed -n ${SLURM_ARRAY_TASK_ID}p`
chr=`awk -F "\t" '{print $3}' $guide_file | sed -n ${SLURM_ARRAY_TASK_ID}p`

#
echo "now processing" ${sample} "at" ${chr} "w/ file==>" ${bamfile}
#

if [[ -e "${bamfile}" ]]
then
	status="exist"
else 
	status="WARNING"
fi

echo $sample $status >> ./${PIPELINE}.completion.log

### extract simplified vcf file per chromosome per sample
  echo "Getting simple vcf"

    bcftools view \
    -s ${sample} \
    -r ${chr} \
    -O v \
    $MASTERVCF \
	-o $WORK_FOLDER/${sample}_${chr}.vcf

	### extract reads from bam file per chromosome per sample
  echo "getting bam"

    samtools view \
    -b \
    -q 20 \
    -M \
    -O BAM \
    -@ $CPU \
    ${bamfile} \
    ${chr} \
     > $WORK_FOLDER/${sample}_${chr}.bam

	## add SM tag and index chr bam file	
	  echo "adding metadata"

	Group_library="SingleInds"
	Library_Platform="illumina"
	Group_platform="SIDMEL"
	
	java -jar $PICARD AddOrReplaceReadGroups \
	I=$WORK_FOLDER/${sample}_${chr}.bam \
	O=$WORK_FOLDER/${sample}_${chr}.SM.bam \
	RGLB=$Group_library \
	RGPL=$Library_Platform \
	RGPU=$Group_platform \
	RGSM=${sample}

	samtools \
	index \
	$WORK_FOLDER/${sample}_${chr}.SM.bam


	### run whatshap
  echo "whatshapp"

  whatshap \
  phase \
  --reference $reference \
  -o $PHASED_FOLDER/${sample}.${chr}.phase.vcf \
  --chromosome ${chr} \
  --sample ${sample} \
  $WORK_FOLDER/${sample}_${chr}.vcf \
  $WORK_FOLDER/${sample}_${chr}.SM.bam

	### clean up
  rm $WORK_FOLDER/${sample}_${chr}.vcf
  rm $WORK_FOLDER/${sample}_${chr}.bam
  rm $WORK_FOLDER/${sample}_${chr}.SM.bam
  rm $WORK_FOLDER/${sample}_${chr}.SM.bam.bai

echo "done" `date`

