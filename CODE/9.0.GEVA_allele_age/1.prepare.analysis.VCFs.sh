#! /bin/bash

#SBATCH -J prep.geva 
#SBATCH --ntasks-per-node=10 # 10 cores
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 
#SBATCH --mem 40G
#SBATCH -p standard
#SBATCH --account jcbnunez

### Prepare VCF files for GEVA ANALYSIS
module load vcftools
module load tabix

#### data in:
input_vcf=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz

############ slice CM
vcftools \
--gzvcf $input_vcf \
--recode \
--recode-INFO-all \
--out CM_pops.2L.all.phased \
--keep ALL_HOMOKAR_CM_samps_OnlyNames.txt \
--chr 2L 

bgzip CM_pops.2L.all.phased.recode.vcf
tabix CM_pops.2L.all.phased.recode.vcf.gz

############ slice CM
vcftools \
--gzvcf $input_vcf \
--recode \
--recode-INFO-all \
--out CM_pops.2L.het.phased \
--keep Het_samps_OnlyNames.txt \
--chr 2L 

bgzip CM_pops.2L.het.phased.recode.vcf
tabix CM_pops.2L.het.phased.recode.vcf.gz


############ slice CM
vcftools \
--gzvcf $input_vcf \
--recode \
--recode-INFO-all \
--out CM_pops.2L.inv.phased \
--keep Inv_samps_OnlyNames.txt \
--chr 2L 

bgzip CM_pops.2L.inv.phased.recode.vcf
tabix CM_pops.2L.inv.phased.recode.vcf.gz


############ slice CM
vcftools \
--gzvcf $input_vcf \
--recode \
--recode-INFO-all \
--out CM_pops.2L.std.phased \
--keep Std_samps_OnlyNames.txt \
--chr 2L 

bgzip CM_pops.2L.std.phased.recode.vcf
tabix CM_pops.2L.std.phased.recode.vcf.gz

date
echo "done"