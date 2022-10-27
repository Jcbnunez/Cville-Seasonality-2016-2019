module load vcftools
module load tabix
module load bcftools

cd  /scratch/yey2sn/Overwintering_ms/7.LD/

input=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz
positions=./matched_controls_SNPids.txt

#bcftools query -f 'pos=%ID\n' $input | head -3

vcftools --gzvcf $input \
--recode --recode-INFO-all \
--out VA_ch.LD.matched_controls.2l \
--snps $positions

bgzip VA_ch.LD.matched_controls.2l.recode.vcf
tabix VA_ch.LD.matched_controls.2l.recode.vcf.gz

#Make the iterator object
bcftools query -f '%ID\n' VA_ch.LD.matched_controls.2l.recode.vcf.gz > ControlMatched.LD.iterator.txt