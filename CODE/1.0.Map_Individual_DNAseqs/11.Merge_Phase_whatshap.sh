#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=170G
#SBATCH --time=72:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez

###########################################################################
#Parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
###########################################################################

### load modules
  module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822

### where are the intermediary phased files
intervals=/scratch/yey2sn/phasing_droso_data/Whatshap/Intervals_for_phasing.txt

INTERMEDIARY_FILES=/scratch/yey2sn/phasing_droso_data/Whatshap/pahsed_output

MERGED_WHATSHAP=/scratch/yey2sn/phasing_droso_data/Whatshap/chr_phased

WD=/scratch/yey2sn/phasing_droso_data/Whatshap/

## make sure we are in the WD
cd $WD

###########################################################################
# Select chromosome
chr=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`

###########################################################################
# Generate folders

if [ -d "MERGED_WHATSHAP" ]
then
	echo "Working mapping_stats folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir MERGED_WHATSHAP
	date
fi

###########################################################################

# bgzip vcf files
for f in $INTERMEDIARY_FILES/*.${chr}.phase.vcf; 
do
  echo "bgzipping File -> $f"
  bgzip \
  -c \
  -@ $CPU \
  -i \
  ${f} > ${f}.gz
done

# index bgzippped files
for f in $INTERMEDIARY_FILES/*.${chr}.phase.vcf.gz; 
do
  echo "indexing File -> $f"
  tabix \
  -p vcf \
  -f \
  ${f}
done

### make file list
  ls -d $INTERMEDIARY_FILES/*.${chr}.phase.vcf.gz \
  > $INTERMEDIARY_FILES/${chr}.list

### Merge files
bcftools \
merge \
-l $INTERMEDIARY_FILES/${chr}.list \
-o  $MERGED_WHATSHAP/${chr}.whatshapp.bcf \
-O b \
--threads $CPU

# index
  bcftools index \
  --threads $CPU \
  $MERGED_WHATSHAP/${chr}.whatshapp.bcf

echo "done" `date`