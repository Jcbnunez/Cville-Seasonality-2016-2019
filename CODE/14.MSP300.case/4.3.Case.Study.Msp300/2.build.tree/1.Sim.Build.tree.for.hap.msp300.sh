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

REF=../dsim-mod.fasta
core.variant.d.sim=4996892
##### 5364427-5364485
start_win=4995892
end_win=4997892
win_name="Mps300.Exon.4995892.4997892"
POS="2L:4995892-4997892"
SP="sim"
### ---> chr2L:4993063-5000175 in D.sim v2
### ---> target 4996892

samtools faidx $REF Dsim_Scf_2L:4995892-4997892 > dsim.2L:4995892-4997892.fasta

## slice 
vcftools \
--gzvcf ../Dsim_all_2L.maf0.01.recode.vcf.gz \
--recode \
--recode-INFO-all \
--out $SP.$win_name \
--chr 2L \
--keep Dsim.samps.for.tree.txt \
--from-bp $start_win \
--to-bp $end_win 

bgzip $SP.$win_name.recode.vcf
tabix $SP.$win_name.recode.vcf.gz

###
mkdir Dsim.haps.f

### neeed to change dsim.2L:4995892-4997892.fasta header to 2L
while read SAMP;
#SAMP=Sz111
do
bcftools consensus \
-H 1 \
-f dsim.2L:4995892-4997892.fasta \
-o Dsim.haps.f/$SAMP.haplotype0.fasta \
-s $SAMP \
$SP.$win_name.recode.vcf.gz

sed -i 's/^>'${POS}'/>'${SAMP}'.0/g' Dsim.haps.f/$SAMP.haplotype0.fasta
done < Dsim.samps.for.tree.txt


