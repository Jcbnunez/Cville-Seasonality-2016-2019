module load vcftools
module load tabix
module load bcftools

input=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz
positions=./VA_ch_0.05_ith0_variant_id.txt

#bcftools query -f 'pos=%ID\n' $input | head -3

vcftools --gzvcf $input \
--recode --recode-INFO-all \
--out VA_ch_0.05_ith0.forLD.2l \
--snps $positions

bgzip VA_ch_0.05_ith0.forLD.2l.recode.vcf
tabix VA_ch_0.05_ith0.forLD.2l.recode.vcf.gz

#Make the iterator object
bcftools query -f '%ID\n' VA_ch_0.05_ith0.forLD.2l.recode.vcf.gz > VA_ch_0.05_ith0.forLD.2l.iterator.txt