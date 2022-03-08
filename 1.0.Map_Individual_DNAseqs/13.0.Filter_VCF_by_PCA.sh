#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --partition=largemem
#SBATCH --account=jcbnunez

###########################################################################
#Parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
###########################################################################

### load modules
module load vcftools
module load tabix
module load bcftools

### where are the intermediary phased files
MERGED_WHATSHAP=/scratch/yey2sn/phasing_droso_data/Whatshap/MERGED_WHATSHAP_noSRs

FILTER=mis80_minQ100_minDP10_maf5_thin1k
PIPELINE=PrePCA
chr=3L

WD=/scratch/yey2sn/phasing_droso_data/Whatshap/

## make sure we are in the WD
cd $WD

###########################################################################
###########################################################################
# Generate folders

if [ -d "PCA_filter" ]
then
	echo "Working mapping_stats folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir PCA_filter
	date
fi

###########################################################################
###########################################################################
# Filter preeliminaru VCF
###########################################################################
###########################################################################

# Filter
vcftools \
--gzvcf $MERGED_WHATSHAP/${chr}.whatshapp.noSR.vcf.gz \
--remove-filtered-all \
--max-missing 0.80 \
--minQ 100 \
--minDP 10 \
--maf 0.05 \
--thin 1000 \
--recode \
--recode-INFO-all \
--out ./PCA_filter/$PIPELINE.$FILTER.${chr}

#bgzip and tabix
bgzip ./PCA_filter/$PIPELINE.$FILTER.${chr}.recode.vcf
tabix ./PCA_filter/$PIPELINE.$FILTER.${chr}.recode.vcf.gz

#Get samps names
bcftools query \
./PCA_filter/$PIPELINE.$FILTER.${chr}.recode.vcf.gz -l snp.vcf \
> ./PCA_filter/$PIPELINE.$FILTER.${chr}.sample_names.txt

echo "done" $(date)
