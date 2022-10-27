#!/usr/bin/env bash
#
#SBATCH -J make_ld_vcf # A single job name for the array
#SBATCH -c 1 
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 4 hours
#SBATCH --mem 4G
#SBATCH -o ./slurmOut/make_ld_vcf.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/make_ld_vcf.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

## filter VCF for informative SNPs
## 

#Load Modules
module load vcftools
module load bcftools
module load tabix
module load picard

####
IN_GZVCF=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz
FILTER=ForPlink
POP=CM

vcftools \
--gzvcf $IN_GZVCF \
--chr 2L \
--maf 0.01 \
--max-alleles 2 \
--recode \
--recode-INFO-all \
--out $FILTER.$POP.2L

bgzip $FILTER.$POP.2L.recode.vcf
tabix $FILTER.$POP.2L.recode.vcf.gz

