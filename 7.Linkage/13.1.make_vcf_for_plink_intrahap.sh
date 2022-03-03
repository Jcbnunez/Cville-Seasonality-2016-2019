module load vcftools
module load tabix
module load bcftools

input=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz
positions=./VA_ch_0.05_ith0_variant_id.txt
inv_inds=./Inv_samps_OnlyNames.txt
#std_invs=


#bcftools query -f 'pos=%ID\n' $input | head -3

vcftools --gzvcf $input \
--recode --recode-INFO-all \
--out VA_ch_0.05_ith0.forLD.2l.INV \
--snps $positions \
--keep $inv_inds \
--max-alleles 2 \
--min-alleles 2 

bgzip VA_ch_0.05_ith0.forLD.2l.INV.recode.vcf
tabix VA_ch_0.05_ith0.forLD.2l.INV.recode.vcf.gz

#Make the iterator object
bcftools query -f '%ID\n' VA_ch_0.05_ith0.forLD.2l.INV.recode.vcf.gz \
 > VA_ch_0.05_ith0.forLD.2l.INV.iterator.txt