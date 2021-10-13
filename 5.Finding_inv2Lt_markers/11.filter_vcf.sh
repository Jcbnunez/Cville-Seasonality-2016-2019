## filter VCF for informative SNPs
## 


#Load Modules
module load vcftools
module load bcftools
module load tabix


IN_GZVCF=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz
  
inf_markers=/project/berglandlab/Dmel_genomic_resources/Inversions/CM_2LT_markers/in2lt_ld_informative_markers.txt

FILTER=inv2Lt_markers

vcftools \
--gzvcf $IN_GZVCF \
--positions $inf_markers \
--recode \
--recode-INFO-all \
--out $FILTER.2L.vcf

#bgzip and tabix
bgzip $PIPELINE.$FILTER.recode.vcf
tabix $PIPELINE.$FILTER.recode.vcf.gz