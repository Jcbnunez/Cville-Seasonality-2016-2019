#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=120G
#SBATCH --time=72:00:00
#SBATCH --partition=largemem
#SBATCH --account=jcbnunez

###########################################################################
#Parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
###########################################################################

### load modules
module load gcc/9.2.0
module load bedtools/2.29.2
module load tabix
module load bcftools

### where are the intermediary phased files
intervals=/scratch/yey2sn/phasing_droso_data/Whatshap/Intervals_for_phasing.txt

MERGED_WHATSHAP=/scratch/yey2sn/phasing_droso_data/Whatshap/MERGED_WHATSHAP

repetitive_regions=/scratch/yey2sn/phasing_droso_data/Whatshap/Repeats.bed

WD=/scratch/yey2sn/phasing_droso_data/Whatshap/

## make sure we are in the WD
cd $WD

###########################################################################
# Select chromosome
chr=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`

###########################################################################
###########################################################################
# Generate folders

if [ -d "MERGED_WHATSHAP_noSRs" ]
then
	echo "Working mapping_stats folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir MERGED_WHATSHAP_noSRs
	date
fi

###########################################################################

## Convert sample to VCF
bcftools view $MERGED_WHATSHAP/${chr}.whatshapp.bcf \
	-O z \
	--threads $CPU \
	-o ./MERGED_WHATSHAP_noSRs/${chr}.whatshapp.vcf \
	
## Run pipeline
  bedtools intersect -v -header \
	 -a ./MERGED_WHATSHAP_noSRs/${chr}.whatshapp.vcf \
     -b <(sort -k1,1 -k2,2n $repetitive_regions) \
	 > ./MERGED_WHATSHAP_noSRs/${chr}.whatshapp.noSR.vcf
	 
## Compress and tabix
	bgzip ./MERGED_WHATSHAP_noSRs/${chr}.whatshapp.noSR.vcf
	tabix ./MERGED_WHATSHAP_noSRs/${chr}.whatshapp.noSR.vcf.gz

## Clean up
rm ./MERGED_WHATSHAP_noSRs/${chr}.whatshapp.vcf

