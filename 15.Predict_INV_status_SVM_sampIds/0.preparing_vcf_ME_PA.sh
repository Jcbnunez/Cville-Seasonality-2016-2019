#!/usr/bin/env bash
#
#SBATCH -J prepVCF # A single job name for the array
#SBATCH --ntasks-per-node=5 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 12:40:00 ### 20 minutes per job per 10 permutations
#SBATCH --mem 25G
#SBATCH -p standard
#SBATCH --account berglandlab_standard


### Extract samples of interest
### 
## filter to only individuals of interest

module load vcftools
module load tabix

SAMPLED_SNPS=/project/berglandlab/Dmel_genomic_resources/Inversions/in2lt_ld_47snps_informative_markers.txt

#taylor specific
Taylor_samps=/scratch/yey2sn/Overwintering_ms/15.ME_PA_extra_pops/nonDGRP_lineids.taylor.txt


###taylor
vcftools \
--gzvcf /project/berglandlab/Dmel_Single_Individuals/Taylors_Data_bams_vcf/Dmel_inds_Taylor.wSNPids.vcf.gz \
--keep $Taylor_samps \
--snps $SAMPLED_SNPS \
--recode \
--recode-INFO-all \
--out Taylors.2L.Rand_INV_Markers

## Compress and tabix
	bgzip Taylors.2L.Rand_INV_Markers.recode.vcf
	tabix Taylors.2L.Rand_INV_Markers.recode.vcf.gz
