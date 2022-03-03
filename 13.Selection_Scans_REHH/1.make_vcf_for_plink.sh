module load vcftools
module load tabix
module load bcftools

#### Make guide object
cat \
/scratch/yey2sn/Overwintering_ms/11.Haplotypes/Inv_samps_OnlyNames.txt \
/scratch/yey2sn/Overwintering_ms/11.Haplotypes/Std_samps_OnlyNames.txt \
> Homozyg_kars.samps.txt


####
input=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz

#bcftools query -f 'pos=%ID\n' $input | head -3

vcftools --gzvcf $input \
--recode --recode-INFO-all \
--chr 2L \
--out VA_ch.HomoKarGenos \
--keep  Homozyg_kars.samps.txt

bgzip VA_ch.HomoKarGenos.recode.vcf
tabix VA_ch.HomoKarGenos.recode.vcf.gz


# Make std homozygous

vcftools --gzvcf $input \
--recode --recode-INFO-all \
--chr 2L \
--out VA_ch.STD.hom \
--keep /scratch/yey2sn/Overwintering_ms/11.Haplotypes/Std_samps_OnlyNames.txt

bgzip VA_ch.STD.hom.recode.vcf
tabix VA_ch.STD.hom.recode.vcf.gz


# Make inv homozygous

vcftools --gzvcf $input \
--recode --recode-INFO-all \
--chr 2L \
--out VA_ch.INV.hom \
--keep /scratch/yey2sn/Overwintering_ms/11.Haplotypes/Inv_samps_OnlyNames.txt

bgzip VA_ch.INV.hom.recode.vcf
tabix VA_ch.INV.hom.recode.vcf.gz

