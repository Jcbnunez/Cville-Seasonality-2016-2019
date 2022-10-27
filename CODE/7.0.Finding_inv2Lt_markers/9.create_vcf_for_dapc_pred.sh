## filter VCF for informative SNPs
## 

#Load Modules
module load vcftools
module load bcftools
module load tabix
module load picard


#bcftools query -f '%CHROM %POS %ID\n' $IN_GZVCF | head -3

inf_markers=/project/berglandlab/Dmel_genomic_resources/Inversions/CM_2LT_markers/in2lt_ld_informative_markers.txt
FILTER=inv2Lt_markers


####
IN_GZVCF=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz
POP=DGRP

vcftools \
--gzvcf $IN_GZVCF \
--snps $inf_markers \
--recode \
--recode-INFO-all \
--out $FILTER.$POP.2L

bgzip $FILTER.$POP.2L.recode.vcf
tabix $FILTER.$POP.2L.recode.vcf.gz

####

####
IN_GZVCF=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz
POP=CM

vcftools \
--gzvcf $IN_GZVCF \
--snps $inf_markers \
--recode \
--recode-INFO-all \
--out $FILTER.$POP.2L

bgzip $FILTER.$POP.2L.recode.vcf
tabix $FILTER.$POP.2L.recode.vcf.gz

