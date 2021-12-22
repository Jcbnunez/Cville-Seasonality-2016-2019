module load vcftools
module load tabix

input=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz
positions=/scratch/yey2sn/Overwintering_ms/7.Haploreconstruct/SNPS_to_select.txt

vcftools --gzvcf $input \
--recode --recode-INFO-all \
--out vcf_for_ld_analysis.2l.vcf \
--positions $positions

bgzip vcf_for_ld_analysis.2l.vcf.recode.vcf
tabix vcf_for_ld_analysis.2l.vcf.recode.vcf.gz