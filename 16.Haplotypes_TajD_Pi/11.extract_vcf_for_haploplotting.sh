module load vcftools
module load tabix
module load bcftools

#### Make guide object
positions=Loci_for_haplotype_plotting_windows.txt

####
cm_vcf2l=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz
dgrp_vcf2l=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz
wolrd_vcf=/project/berglandlab/Dmel_Single_Individuals/Taylors_Data_bams_vcf/Dmel_inds_Taylor.wSNPids.vcf.gz


#bcftools query -f 'pos=%ID\n' $input | head -3

###
cat \
/scratch/yey2sn/Overwintering_ms/16.Haplotypes/CM.INV.samps.txt \
/scratch/yey2sn/Overwintering_ms/16.Haplotypes/CM.STD.samps.txt \
> CM.Homozyg_kars.samps.txt

cat \
/scratch/yey2sn/Overwintering_ms/16.Haplotypes/INV_DGRP_OnlyNames.txt \
/scratch/yey2sn/Overwintering_ms/16.Haplotypes/STD_DGRP_OnlyNames.txt \
> DGRP.Homozyg_kars.samps.txt

########
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

########
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

########
vcftools --gzvcf $wolrd_vcf \
--recode --recode-INFO-all \
--chr 2L \
--out Taylor.haplo.snps \
--max-alleles  2 \
--remove-indels \
--snps $positions \
--keep /project/berglandlab/Dmel_Single_Individuals/in2Lt_status_info/samples_hom_In2lt.ids.Taylor.txt

bgzip Taylor.haplo.snps.recode.vcf
tabix Taylor.haplo.snps.recode.vcf.gz

#######
#######
#######
#######
bcftools merge \
CM.haplo.snps.recode.vcf.gz \
DGRP.haplo.snps.recode.vcf.gz \
Taylor.haplo.snps.recode.vcf.gz \
> MC.DGRP.Taylor.merged.vcf

bgzip MC.DGRP.Taylor.merged.vcf
tabix MC.DGRP.Taylor.merged.vcf.gz
#######
########
#######
########


# filter ready for import 
vcftools --gzvcf MC.DGRP.Taylor.merged.vcf.gz \
--recode --recode-INFO-all \
--chr 2L \
--out MC.DGRP.Taylor.merged.readyForImport \
--max-alleles  2 \
--remove-indels 
bgzip MC.DGRP.Taylor.merged.readyForImport.recode.vcf
tabix MC.DGRP.Taylor.merged.readyForImport.recode.vcf.gz

### filter maf for haplonet
#=# vcftools --gzvcf MC.DGRP.merged.readyForImport.recode.vcf.gz \
#=# --recode --recode-INFO-all \
#=# --chr 2L \
#=# --out MC.DGRP.merged.forHap.maf \
#=# --max-alleles  2 \
#=# --maf 0.05 \
#=# --remove-indels 
#=# bgzip MC.DGRP.merged.forHap.maf.recode.vcf
#=# tabix MC.DGRP.merged.forHap.maf.recode.vcf.gz

##
##
##
##
##
##
## Process the heterozygous

#=# vcftools --gzvcf $cm_vcf2l \
#=# --recode --recode-INFO-all \
#=# --chr 2L \
#=# --out CM.het.snps \
#=# --max-alleles  2 \
#=# --snps $positions \
#=# --remove CM.Homozyg_kars.samps.txt
#=# 
#=# bgzip CM.het.snps.recode.vcf
#=# tabix CM.het.snps.recode.vcf.gz
#=# 
#=# vcftools --gzvcf $dgrp_vcf2l \
#=# --recode --recode-INFO-all \
#=# --chr 2L \
#=# --out DGRP.het.snps \
#=# --max-alleles  2 \
#=# --snps $positions \
#=# --remove DGRP.Homozyg_kars.samps.txt
#=# 
#=# bgzip DGRP.het.snps.recode.vcf
#=# tabix DGRP.het.snps.recode.vcf.gz
#=# 
#=# bcftools merge CM.het.snps.recode.vcf.gz DGRP.het.snps.recode.vcf.gz > MC.DGRP.het.merged.vcf
#=# bgzip MC.DGRP.het.merged.vcf
#=# tabix MC.DGRP.het.merged.vcf.gz
#=# 
#=# # filter ready for import 
#=# vcftools --gzvcf MC.DGRP.het.merged.vcf.gz \
#=# --recode --recode-INFO-all \
#=# --chr 2L \
#=# --out MC.DGRP.het.merged.readyForImport \
#=# --max-alleles  2 \
#=# --remove-indels 
#=# bgzip MC.DGRP.het.merged.readyForImport.recode.vcf
#=# tabix MC.DGRP.het.merged.readyForImport.recode.vcf.gz




