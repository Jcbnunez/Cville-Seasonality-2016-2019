module load vcftools
module load samtools
module load tabix
module load bcftools
module load picard

#load in the reference genome
reference=/project/berglandlab/Dmel_genomic_resources/References/holo_dmel_6.12.fa

# load in the gavcf
input_vcf=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz
dgrp_vcf2l=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz

start_win=6.00e6
end_win=6.75e6
win_name="win6"
POS="2L:6000000-6750000"

## slice CM
vcftools \
--gzvcf $input_vcf \
--recode \
--recode-INFO-all \
--out $win_name.CM_pops.2L.invReg.phased \
--chr 2L 
--from-bp $start_win \
--to-bp $end_win 

bgzip CM_pops.2L.invReg.phased.recode.vcf
tabix CM_pops.2L.invReg.phased.recode.vcf.gz

## slice DGRP
vcftools \
--gzvcf $dgrp_vcf2l \
--recode \
--recode-INFO-all \
--out $win_name.CM_pops.2L.invReg.phased \
--chr 2L \
--from-bp $start_win \
--to-bp $end_win 

bgzip  $win_name.DGRP.2L.invReg.phased.recode.vcf
tabix  $win_name.DGRP.2L.invReg.phased.recode.vcf.gz


##Merge the 2l vcfs
vcf-merge CM_pops.2L.phased.recode.vcf.gz $dgrp_vcf2l > CM.DGRP.2L.merged.vcf


### guides
retain_loci=Retain_loci_posOnly.txt
retain_samps=Homozyg_samps_OnlyNames.txt

#make ref

samtools faidx $reference $POS >	locus.$POS.2Lt.fasta

#filter the VCF
vcftools \
--gzvcf $input_vcf \
--recode \
--recode-INFO-all \
--out CM_inv2lt_GLM_LD \
#--positions $retain_loci #\
#--keep $retain_samps
#--chr 2L \
#--from-bp 2051600 \
#--to-bp 13204680 \

bgzip CM_inv2lt_GLM_LD.recode.vcf
tabix CM_inv2lt_GLM_LD.recode.vcf.gz


#filter the VCF
vcftools \
--gzvcf $dgrp_vcf2l \
--recode \
--recode-INFO-all \
--out DGRP_inv2lt_GLM_LD \
--positions $retain_loci #\
--keep $retain_samps
--chr 2L \
#--from-bp 2051600 \
#--to-bp 13204680 \

bgzip DGRP_inv2lt_GLM_LD.recode.vcf
tabix DGRP_inv2lt_GLM_LD.recode.vcf.gz



#
#
#
#
#
#
#
#
#
#
#
#
##########
#
#### process homozygotes
#
#
#mkdir Stds_haps
#
#while read SAMP;
#do
#
#echo "now processing" $SAMP
#
#echo "extract 1"
#bcftools consensus \
#-H 1 \
#-f locus.$POS.2Lt.fasta \
#-o $SAMP.STD.GLM_LD.haplotype0.fasta \
#-s $SAMP \
#HomozygSamps_GLM_LD.recode.vcf.gz
#
#echo "extract 2"
#bcftools consensus \
#-H 2 \
#-f locus.$POS.2Lt.fasta \
#-o $SAMP.STD.GLM_LD.haplotype1.fasta \
#-s $SAMP \
#HomozygSamps_GLM_LD.recode.vcf.gz
#
#echo "parse 2"
#sed -e 's/^>'${POS}'/>'${SAMP}'.0/g' $SAMP.STD.GLM_LD.haplotype0.fasta > $SAMP.STD.GLM_LD.haplotype0.named.fasta
#rm $SAMP.STD.GLM_LD.haplotype0.fasta
#mv $SAMP.STD.GLM_LD.haplotype0.named.fasta ./Stds_haps
#
#echo "parse 2"
#sed -e 's/^>'${POS}'/>'${SAMP}'.1/g' $SAMP.STD.GLM_LD.haplotype1.fasta > $SAMP.STD.GLM_LD.haplotype1.named.fasta
#rm $SAMP.STD.GLM_LD.haplotype1.fasta
#mv $SAMP.STD.GLM_LD.haplotype1.named.fasta ./Stds_haps
#
#done < $retain_samps
#
#cat ./Stds_haps/*.named.fasta > STD_INV.CM.joint.fasta
#
#################
#################
#################
#################
#### Extract hets
#retain_hets=Hetg_samps_OnlyNames.txt
#
##filter the VCF
#vcftools \
#--gzvcf $input_vcf \
#--recode \
#--recode-INFO-all \
#--out HetSamps_GLM_LD \
#--positions $retain_loci \
#--keep $retain_hets
#
#bgzip HetSamps_GLM_LD.recode.vcf
#tabix HetSamps_GLM_LD.recode.vcf.gz
#
### process hets
#mkdir het_haps
#
#while read SAMP;
#do
#
#echo "now processing" $SAMP
#
#echo "extract 1"
#bcftools consensus \
#-H 1 \
#-f locus.$POS.2Lt.fasta \
#-o $SAMP.HET.GLM_LD.haplotype0.fasta \
#-s $SAMP \
#HetSamps_GLM_LD.recode.vcf.gz
#
#echo "extract 2"
#bcftools consensus \
#-H 2 \
#-f locus.$POS.2Lt.fasta \
#-o $SAMP.HET.GLM_LD.haplotype1.fasta \
#-s $SAMP \
#HetSamps_GLM_LD.recode.vcf.gz
#
#echo "parse 2"
#sed -e 's/^>'${POS}'/>'${SAMP}'.0/g' $SAMP.HET.GLM_LD.haplotype0.fasta > $SAMP.HET.GLM_LD.haplotype0.named.fasta
#rm $SAMP.HET.GLM_LD.haplotype0.fasta
#mv $SAMP.HET.GLM_LD.haplotype0.named.fasta ./het_haps
#
#echo "parse 2"
#sed -e 's/^>'${POS}'/>'${SAMP}'.1/g' $SAMP.HET.GLM_LD.haplotype1.fasta > $SAMP.HET.GLM_LD.haplotype1.named.fasta
#rm $SAMP.HET.GLM_LD.haplotype1.fasta
#mv $SAMP.HET.GLM_LD.haplotype1.named.fasta ./het_haps
#
#done < $retain_hets
#
#cat ./het_haps/*.named.fasta > HET.CM.joint.fasta
#

