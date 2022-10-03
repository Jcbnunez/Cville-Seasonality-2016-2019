### 7. prepare vcf for visualizing haplotypes

module load vcftools
module load tabix

#load in the reference genome
reference=/project/berglandlab/Dmel_genomic_resources/References/holo_dmel_6.12.fa

# load in the gavcf
input_vcf=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz
dgrp_vcf2l=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz


## slice CM
vcftools \
--gzvcf $input_vcf \
--recode \
--recode-INFO-all \
--out CM_pops.2L.invReg.phased \
--chr 2L 
--from-bp 2000000 \
--to-bp 13210000 

bgzip CM_pops.2L.invReg.phased.recode.vcf
tabix CM_pops.2L.invReg.phased.recode.vcf.gz

## slice DGRP
vcftools \
--gzvcf $dgrp_vcf2l \
--recode \
--recode-INFO-all \
--out DGRP.2L.invReg.phased \
--chr 2L \
--from-bp 2000000 \
--to-bp 13210000 

bgzip DGRP.2L.invReg.phased.recode.vcf
tabix DGRP.2L.invReg.phased.recode.vcf.gz
