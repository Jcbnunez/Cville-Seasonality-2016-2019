module load vcftools
module load tabix
module load bcftools

#### Make guide object
positions=Loci_for_haplotype_plotting_windows.txt

####
cm_vcf2l=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz
dgrp_vcf2l=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz

#bcftools query -f 'pos=%ID\n' $input | head -3

###
cat \
/scratch/yey2sn/Overwintering_ms/11.Haplotypes/Inv_samps_OnlyNames.txt \
/scratch/yey2sn/Overwintering_ms/11.Haplotypes/Std_samps_OnlyNames.txt \
> CM.Homozyg_kars.samps.txt

cat \
/scratch/yey2sn/Overwintering_ms/11.Haplotypes/INV_DGRP_OnlyNames.txt \
/scratch/yey2sn/Overwintering_ms/11.Haplotypes/STD_DGRP_OnlyNames.txt \
> DGRP.Homozyg_kars.samps.txt

##

vcftools --gzvcf $cm_vcf2l \
--recode --recode-INFO-all \
--chr 2L \
--out CM.haplo.snps \
--max-alleles  2 \
--remove-indels \
--snps $positions \
--keep CM.Homozyg_kars.samps.txt

bgzip CM.haplo.snps.recode.vcf
tabix CM.haplo.snps.recode.vcf.gz

vcftools --gzvcf $dgrp_vcf2l \
--recode --recode-INFO-all \
--chr 2L \
--out DGRP.haplo.snps \
--max-alleles  2 \
--remove-indels \
--snps $positions \
--keep DGRP.Homozyg_kars.samps.txt

bgzip DGRP.haplo.snps.recode.vcf
tabix DGRP.haplo.snps.recode.vcf.gz

bcftools merge CM.haplo.snps.recode.vcf.gz DGRP.haplo.snps.recode.vcf.gz > MC.DGRP.merged.vcf
bgzip MC.DGRP.merged.vcf
tabix MC.DGRP.merged.vcf.gz

# filter ready for import 
vcftools --gzvcf MC.DGRP.merged.vcf.gz \
--recode --recode-INFO-all \
--chr 2L \
--out MC.DGRP.merged.readyForImport \
--max-alleles  2 \
--remove-indels 
bgzip MC.DGRP.merged.readyForImport.recode.vcf
tabix MC.DGRP.merged.readyForImport.recode.vcf.gz

### filter maf for haplonet
vcftools --gzvcf MC.DGRP.merged.readyForImport.recode.vcf.gz \
--recode --recode-INFO-all \
--chr 2L \
--out MC.DGRP.merged.forHap.maf \
--max-alleles  2 \
--maf 0.05 \
--remove-indels 
bgzip MC.DGRP.merged.forHap.maf.recode.vcf
tabix MC.DGRP.merged.forHap.maf.recode.vcf.gz

##
##
##
##
##
##
## Process the heterozygous

vcftools --gzvcf $cm_vcf2l \
--recode --recode-INFO-all \
--chr 2L \
--out CM.het.snps \
--max-alleles  2 \
--snps $positions \
--remove CM.Homozyg_kars.samps.txt

bgzip CM.het.snps.recode.vcf
tabix CM.het.snps.recode.vcf.gz

vcftools --gzvcf $dgrp_vcf2l \
--recode --recode-INFO-all \
--chr 2L \
--out DGRP.het.snps \
--max-alleles  2 \
--snps $positions \
--remove DGRP.Homozyg_kars.samps.txt

bgzip DGRP.het.snps.recode.vcf
tabix DGRP.het.snps.recode.vcf.gz

bcftools merge CM.het.snps.recode.vcf.gz DGRP.het.snps.recode.vcf.gz > MC.DGRP.het.merged.vcf
bgzip MC.DGRP.het.merged.vcf
tabix MC.DGRP.het.merged.vcf.gz

# filter ready for import 
vcftools --gzvcf MC.DGRP.het.merged.vcf.gz \
--recode --recode-INFO-all \
--chr 2L \
--out MC.DGRP.het.merged.readyForImport \
--max-alleles  2 \
--remove-indels 
bgzip MC.DGRP.het.merged.readyForImport.recode.vcf
tabix MC.DGRP.het.merged.readyForImport.recode.vcf.gz




