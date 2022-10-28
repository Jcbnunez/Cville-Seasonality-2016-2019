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
start_win=5092177
end_win=5292177
win_name="Mps300.betascan"
POS="2L:5092177-5292177"
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
--from-bp $start_win \
--to-bp $end_win 

bgzip $SP.$win_name.recode.vcf
tabix $SP.$win_name.recode.vcf.gz

bcftools query -l $SP.$win_name.recode.vcf.gz > $SP.$win_name.samp.names.txt

###
mkdir Dmel.haps.betascan0

while read SAMP;
do
echo $SAMP;
bcftools consensus \
-H 1 \
-f $SP.$POS.fasta \
-o  Dmel.haps.betascan0/$SAMP.haplotype0.fasta \
-s $SAMP \
$SP.$win_name.recode.vcf.gz

sed -i 's/^>'${POS}'/>'${SP}.${SAMP}'.0/g'  Dmel.haps.betascan0/$SAMP.haplotype0.fasta
done < $SP.$win_name.samp.names.txt

###
mkdir Dmel.haps.betascan1

while read SAMP;
do
echo $SAMP;
bcftools consensus \
-H 1 \
-f $SP.$POS.fasta \
-o  Dmel.haps.betascan1/$SAMP.haplotype0.fasta \
-s $SAMP \
$SP.$win_name.recode.vcf.gz

sed -i 's/^>'${POS}'/>'${SP}.${SAMP}'.1/g'  Dmel.haps.betascan1/$SAMP.haplotype0.fasta
done < $SP.$win_name.samp.names.txt

