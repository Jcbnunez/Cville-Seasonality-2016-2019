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

REF=/project/berglandlab/Dmel_genomic_resources/Other_drosophilids/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_genomic.fna
core.variant.d.sim=9456177
##### 5364427-5364485
start_win=9455179
end_win=9457180
win_name="Mps300.Exon.9455179.9457180"
POS='NC_052527.2:9455179-9457180'
SP="yak"
### ---> chr2L:4993063-5000175 in D.sim v2
### ---> target 4996892


samtools faidx $REF
samtools faidx $REF $POS > $SP.$POS.fasta


vcf=/scratch/yey2sn/Overwintering_ms/msp300.case/out.vcfs/NC_052527.2.genotyped.raw.vcf.gz
## slice 
vcftools \
--gzvcf $vcf \
--recode \
--recode-INFO-all \
--out $SP.$win_name \
--chr NC_052527.2 \
--keep Dyak.samps.for.tree.txt \
--from-bp $start_win \
--to-bp $end_win 

bgzip $SP.$win_name.recode.vcf
tabix $SP.$win_name.recode.vcf.gz

###
mkdir Dyak.haps.f

while read SAMP;
#SAMP=Sz111
do
bcftools consensus \
-H 1 \
-f $SP.$POS.fasta \
-o Dyak.haps.f/$SAMP.haplotype0.fasta \
-s $SAMP \
$SP.$win_name.recode.vcf.gz

sed -i 's/^>'${POS}'/>'${SAMP}'.0/g' Dyak.haps.f/$SAMP.haplotype0.fasta
done <  Dyak.samps.for.tree.txt


