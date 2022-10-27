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


wd=/scratch/yey2sn/Overwintering_ms/Fig6_analyze_temp_snps
GVZF=ForPlink.CM.2L.recode

# add annotation
bcftools query \
-f '%CHROM\t%POS\t%POS\t%CHROM\_%POS\_SNP\n' \
$wd/$GVZF.vcf.gz > $GVZF.annotation.txt

#Index the annotation
bgzip $GVZF.annotation.txt
tabix -s1 -b2 -e2 $GVZF.annotation.txt.gz

# Add annotation
bcftools annotate \
-a  $GVZF.annotation.txt.gz \
-c CHROM,FROM,TO,ID \
$wd/$GVZF.vcf.gz > $wd/$GVZF.wSNPids.vcf

#tabix
bgzip $wd/$GVZF.wSNPids.vcf
tabix $wd/$GVZF.wSNPids.vcf.gz

#sanity check
bcftools query \
-f '%CHROM\t%POS\t%ID\n' \
$wd/$GVZF.wSNPids.vcf.gz | head
