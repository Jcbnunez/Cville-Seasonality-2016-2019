#### Build tree for Msp300 -- local tree
#### 

#####
#####

module load vcftools
module load samtools
module load tabix
module load bcftools
module load picard
module load tabix
module load bcftools

### bgzipping D sim genome
#bgzip simulans_multisamp_all_chr.vcf
#tabix simulans_multisamp_all_chr.vcf.gz

REF=/project/berglandlab/Dmel_genomic_resources/References/holo_dmel_6.12.fa
core.variant.d.sim=5192177
##### 5364427-5364485
start_win=5191177
end_win=5193177
win_name="Mps300.Exon.5191177.5193177"
POS="2L:5191177-5193177"
SP="mel"

samtools faidx $REF $POS > $SP.$POS.fasta

input_vcf=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz

## slice 
vcftools \
--gzvcf $input_vcf \
--recode \
--recode-INFO-all \
--out $SP.$win_name \
--chr 2L \
--keep Dmel.samps.for.tree.txt \
--from-bp $start_win \
--to-bp $end_win 

bgzip $SP.$win_name.recode.vcf
tabix $SP.$win_name.recode.vcf.gz

###
mkdir Dmel.haps.f

while read SAMP;
#SAMP=Sz111
do
bcftools consensus \
-H 1 \
-f $SP.$POS.fasta \
-o Dmel.haps.f/$SAMP.haplotype0.fasta \
-s $SAMP \
$SP.$win_name.recode.vcf.gz

sed -i 's/^>'${POS}'/>'${SAMP}'.0/g' Dmel.haps.f/$SAMP.haplotype0.fasta
done < Dmel.samps.for.tree.txt


