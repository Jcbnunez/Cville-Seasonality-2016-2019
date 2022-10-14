#!/bin/sh
#
#SBATCH -J add_annot_snp
#SBATCH -c 10
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00 
#SBATCH --mem 40G
#SBATCH -p standard
#SBATCH --account jcbnunez

#load modules
module load tabix
module load bcftools
module load vcftools


vcf=DPGP3.2L.chrfix.recode.vcf

# add annotation
bcftools query \
-f '%CHROM\t%POS\t%POS\t%CHROM\_%POS\_SNP\n' \
$vcf > dpgp3.2L.annotation.txt

#Index the annotation
bgzip dpgp3.2L.annotation.txt
tabix -s1 -b2 -e2 dpgp3.2L.annotation.txt.gz

# Add annotation
bcftools annotate \
-a  dpgp3.2L.annotation.txt.gz \
-c CHROM,FROM,TO,ID \
DPGP3.2L.chrfix.recode.vcf > DPGP3.2L.chrfixAnnot.recode.vcf

#tabix
bgzip DPGP3.2L.chrfixAnnot.recode.vcf
tabix DPGP3.2L.chrfixAnnot.recode.vcf.gz

##
DGRP=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz

bcftools merge \
$DGRP \
DPGP3.2L.chrfixAnnot.recode.vcf.gz \
> DPGP3.DGRP.2L.merged.vcf

bgzip DPGP3.DGRP.2L.merged.vcf
tabix DPGP3.DGRP.2L.merged.vcf.gz

### filter
vcftools \
--gzvcf DPGP3.DGRP.2L.merged.vcf.gz \
--recode --recode-INFO-all \
--max-alleles 2 \
--max-missing 0.8 \
--out DPGP3.DGRP.2L.merged.flt

bgzip DPGP3.DGRP.2L.merged.flt.recode.vcf
tabix DPGP3.DGRP.2L.merged.flt.recode.vcf.gz

