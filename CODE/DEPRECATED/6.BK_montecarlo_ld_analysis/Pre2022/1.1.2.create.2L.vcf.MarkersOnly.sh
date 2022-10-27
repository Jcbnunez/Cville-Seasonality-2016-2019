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
# --> the file "combined_temp_inv_markers.txt" is made by running a cat command between the inversion markers see folder #5 and the temperature from the temperature GLM.. SNPs see folder #4

####
IN_GZVCF=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz
FILTER=ForPlinkMarkersOnly
POP=CM
snp_file=./combined_temp_inv_markers.txt

vcftools \
--gzvcf $IN_GZVCF \
--snps $snp_file \
--recode \
--recode-INFO-all \
--out $FILTER.$POP.2L

bgzip $FILTER.$POP.2L.recode.vcf
tabix $FILTER.$POP.2L.recode.vcf.gz

