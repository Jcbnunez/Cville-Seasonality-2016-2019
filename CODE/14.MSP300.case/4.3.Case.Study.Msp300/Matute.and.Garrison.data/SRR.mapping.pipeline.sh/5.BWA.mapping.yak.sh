#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=8:00:00
#SBATCH -o ./slurmOut/map%A_%a.out # Standard output
#SBATCH -e ./slurmOut/map%A_%a.err # Standard error
#SBATCH --partition=standard
#SBATCH --account=jcbnunez
#SBATCH --array=2-36

####################
#Define parameters
CPU=1
#CPU=$SLURM_CPUS_ON_NODE
#echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
QUAL=40 # Quality threshold for samtools
JAVAMEM=18g # Java memory
####################


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
RAW_READS=/scratch/yey2sn/Overwintering_ms/msp300.case/fasta.files.fld
#Working folder is core folder where this pipeline is being run. -- #out-folder
OUT=/scratch/yey2sn/Overwintering_ms/msp300.case/mapping_output
#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/project/berglandlab/Dmel_genomic_resources/Other_drosophilids/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_genomic.fna

###
Metadat=yak.metadat.mapping.txt

###
sampleId=`sed '1d' $Metadat | awk -F "\t" '{print $1}' | sed -n ${SLURM_ARRAY_TASK_ID}p`
config=`sed '1d' $Metadat | awk -F "\t" '{print $2}' | sed -n ${SLURM_ARRAY_TASK_ID}p`
sp=`sed '1d' $Metadat | awk -F "\t" '{print $3}' | sed -n ${SLURM_ARRAY_TASK_ID}p`
SRR=`sed '1d' $Metadat | awk -F "\t" '{print $4}' | sed -n ${SLURM_ARRAY_TASK_ID}p`
ReadF=`sed '1d' $Metadat | awk -F "\t" '{print $5}' | sed -n ${SLURM_ARRAY_TASK_ID}p`
ReadR=`sed '1d' $Metadat | awk -F "\t" '{print $6}' | sed -n ${SLURM_ARRAY_TASK_ID}p`

echo $sampleId "|" $config "|" $sp "|"  $SRR "|" $ReadF "|" $ReadR "|" 
echo  $CPU "|" $REFERENCE
#ls $RAW_READS
echo $RAW_READS/$ReadF
more $RAW_READS/$ReadF

###
mkdir $OUT/$sp
mkdir $OUT/$sp/$SRR

if [[ $config == 1 ]]; then
 echo "seems this is single end data, lets map it"
 bwa mem \
 -M \
 -t $CPU \
 $REFERENCE \
 $RAW_READS/$ReadF \
 > $OUT/$sp/$SRR/$sampleId.sam

elif [[ $config == 2 ]]; then
 echo "seems this is paired end data, lets map it using a 1-2 approach"
 bwa mem \
 -M \
 -t $CPU \
 $REFERENCE \
 $RAW_READS/$ReadF \
 $RAW_READS/$ReadR \
 > $OUT/$sp/$SRR/$sampleId.sam

else
 echo "I cant tell what type of data is this -- WARNING!"
 echo ${SRR} "Something is wrong at the mapping stage" $(date)
fi

echo "mapping done"
echo "now building bams"

#I will now extract some summary stats
 samtools flagstat \
 --threads $CPU \
 $OUT/$sp/$SRR/$sampleId.sam \
 > $OUT/$sp/$SRR/$sampleId.stats.txt

#build bam files
 samtools view \
 -b \
 -q $QUAL \
 --threads $CPU  \
 $OUT/$sp/$SRR/$sampleId.sam \
 > $OUT/$sp/$SRR/$sampleId.bam

####
####

# Sort with picard
 # Notice that once a file has been sorted it is added the "srt" suffix
 java -Xmx$JAVAMEM \
 -jar $PICARD SortSam \
 I=$OUT/$sp/$SRR/$sampleId.bam \
 O=$OUT/$sp/$SRR/$sampleId.srt.bam \
 SO=coordinate \
 VALIDATION_STRINGENCY=SILENT

#Remove duplicates with picard
 # Notice that once a file has been sorted it is added the "rmdp" suffix
 java -Xmx$JAVAMEM \
 -jar $PICARD MarkDuplicates \
 I=$OUT/$sp/$SRR/$sampleId.srt.bam \
 O=$OUT/$sp/$SRR/$sampleId.srt.rmdp.bam \
 M=$OUT/$sp/$SRR/$sampleId.dupstat.txt \
 VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

#Lets do QC on the bam file
 qualimap bamqc \
 -bam $OUT/$sp/$SRR/$sampleId.srt.rmdp.bam  \
 -outdir $OUT/$sp/$SRR/Qualimap_${SRR} \
 --java-mem-size=$JAVAMEM

## housekeeping
 rm $OUT/$sp/$SRR/$sampleId.sam
 rm $OUT/$sp/$SRR/$sampleId.bam
 rm $OUT/$sp/$SRR/$sampleId.srt.bam

echo "done"
