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

REF=/project/berglandlab/Dmel_genomic_resources/Other_drosophilids/GCF_004382195.2_Dsechelia_ASM438219v2_genomic.fna
core.variant.d.sim=5099990
##### 5364427-5364485
start_win=5098992
end_win=5100990
win_name="Mps300.Exon.5098992.5100990"
POS='NC_045949.1:5098992-5100990'
SP="sech"
### ---> chr2L:4993063-5000175 in D.sim v2
### ---> target 4996892


samtools faidx $REF
samtools faidx $REF $POS > $SP.$POS.fasta


vcf=/scratch/yey2sn/Overwintering_ms/msp300.case/out.vcfs/NC_045949.1.genotyped.raw.vcf.gz
## slice 
vcftools \
--gzvcf $vcf \
--recode \
--recode-INFO-all \
--out $SP.$win_name \
--chr NC_045949.1 \
--keep Dsech.samps.for.tree.txt \
--from-bp $start_win \
--to-bp $end_win 

bgzip $SP.$win_name.recode.vcf
tabix $SP.$win_name.recode.vcf.gz

###
mkdir Dsech.haps.f

while read SAMP;
#SAMP=Sz111
do
bcftools consensus \
-H 1 \
-f $SP.$POS.fasta \
-o Dsech.haps.f/$SAMP.haplotype0.fasta \
-s $SAMP \
$SP.$win_name.recode.vcf.gz

sed -i 's/^>'${POS}'/>'${SAMP}'.0/g' Dsech.haps.f/$SAMP.haplotype0.fasta
done <  Dsech.samps.for.tree.txt


