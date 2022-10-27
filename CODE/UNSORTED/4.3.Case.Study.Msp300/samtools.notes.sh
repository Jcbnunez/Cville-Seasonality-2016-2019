#### Characterize 2L_5192177_SNP in the DGRP.
#### 
#### ---> look for motiff *CCTCCTGGAGAT* -- pos in exon 3835

##### ----> in D. simulans pos is 5156694
#A Large Panel of Drosophila simulans Reveals an Abundance of Common Variants
#Sarah A Signor, Felicia N New, Sergey Nuzhdin Author Notes
#Genome Biology and Evolution, Volume 10, Issue 1, January 2018, Pages 189–206, https://doi.org/10.1093/gbe/evx262
## in sim vs pos should be "4996891"



#Fine scale mapping of genomic introgressions within the Drosophila yakuba clade
#    David A. Turissini,
#    Daniel R. Matute


# ----> get data
wget http://example.com/sample.php https://zenodo.org/record/154261/files/simulans_multisamp_all_chr.vcf.zip
cp /project/berglandlab/courtney/simCline/refgenomes/simulans/dsim-mod.fasta ./


module load exonerate
target=./dsim-mod.fasta
query=./dmel.MSP300.2L.5188342.5195460.AAseq.fasta

exonerate --model protein2genome $query $target > AAmsp300.exo.search.txt

samtools faidx $target Dsim_Scf_2L:4993062-5000175 > dsim.2L.4993062.5000175.msp300.fasta

### Exons of MSP300

module load vcftools
module load samtools
module load tabix
module load bcftools
module load picard
module load tabix
module load bcftools

## Exon of interest in MSP 300 is 2L:5188342-5195460

REF=/project/berglandlab/Dmel_genomic_resources/References/holo_dmel_6.12.fa
samtools faidx $REF 2L:5188342-5195460 > dmel.2L.5188342.5195460.fasta


#Description: fatty acid synthase 1, transcript variant A
samtools faidx $REF 2L:3056598-3068177 > dmel.2L:3056598-3068177.fasta


### Notes
#Drosophila simulans strain w501 chromosome 2L, whole genome shotgun sequence
#Sequence ID of D.sim w501: CM028640.2
#from 5152865 to 5159939
## end of notes

#Description: fatty acid synthase 1, transcript variant A
#CM028640.2:3057155-3068729 
#Drosophila simulans strain w501 chromosome 2L, whole genome shotgun sequence
## 


#load in the reference genome
reference=/project/berglandlab/Dmel_genomic_resources/References/holo_dmel_6.12.fa

# load in the gavcf
input_vcf=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz
dgrp_vcf2l=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz
taylor_vcf=/project/berglandlab/Dmel_Single_Individuals/Taylors_Data_bams_vcf/Dmel_inds_Taylor.wSNPids.vcf.gz

start_win=5188342
end_win=5195460
win_name="Mps300.Exon"
POS="2L:5188342-5195460"

## slice CM
vcftools \
--gzvcf $input_vcf \
--recode \
--recode-INFO-all \
--out $win_name.CM_pops.2L.invReg.phased \
--chr 2L \
--from-bp $start_win \
--to-bp $end_win 
bgzip $win_name.CM_pops.2L.invReg.phased.recode.vcf
tabix $win_name.CM_pops.2L.invReg.phased.recode.vcf.gz

## slice AFRICA
vcftools \
--gzvcf $taylor_vcf \
--recode \
--recode-INFO-all \
--out $win_name.AFR.2L.invReg.phased \
--chr 2L \
--from-bp $start_win \
--to-bp $end_win 
bgzip $win_name.AFR.2L.invReg.phased.recode.vcf
tabix $win_name.AFR.2L.invReg.phased.recode.vcf.gz





####
mkdir msp.case.f
while read SAMP;
do
bcftools consensus \
-H 1 \
-f dmel.2L.5188342.5195460.fasta \
-o msp.case.f/$SAMP.STD.GLM_LD.haplotype0.fasta \
-s $SAMP \
$win_name.CM_pops.2L.invReg.phased.recode.vcf.gz

sed -i 's/^>'${POS}'/>'${SAMP}'.0/g' msp.case.f/$SAMP.STD.GLM_LD.haplotype0.fasta
done < test.samps.homs.txt

while read SAMP;
do
bcftools consensus \
-H 1 \
-f dmel.2L.5188342.5195460.fasta \
-o msp.case.f/$SAMP.STD.GLM_LD.haplotype0.fasta \
-s $SAMP \
$win_name.AFR.2L.invReg.phased.recode.vcf.gz

sed -i 's/^>'${POS}'/>'${SAMP}'.0/g' msp.case.f/$SAMP.STD.GLM_LD.haplotype0.fasta
done < africa.samps.txt


cat msp.case.f/*.fasta \
w501.MSP300.Exon.Dsimulans.fasta \
dyak.MSP300.Exon.Dyak.fasta \
schel.MSP300.Exon.fasta \
> test.haps.fasta

###
module load gcc/9.2.0  openmpi/3.1.6
module load mafft
module load iqtree

mafft test.haps.fasta > test.haps.mafft.fasta

mkdir tree
iqtree -s test.haps.mafft.fasta -bb 1000 -o Dyak.msp300 --prefix tree/msp300


### R
module load gcc/7.1.0
module load openmpi/3.1.4
module load R/4.1.1
module load gdal
module load proj
R

## in R
## 

library("treeio")
library("ggtree")
library("vroom")
library("TDbook")
library("adegenet")
library("FactoMineR")
library("tidyverse")


nwk <- "msp.tree/msp300.treefile"
tree <- read.tree(nwk)
#genotype <- as.data.frame(vroom("tree.metadat.txt"))
#genotype.tip = as.data.frame(genotype[,c("kar")])
#row.names(genotype.tip) = paste(genotype$tip, ".0", sep = "")
dna.seq <- fasta2DNAbin("Msp.300.haps.al.fasta", quiet=FALSE, chunkSize=10, snpOnly=FALSE)

bin.dna <- DNAbin2genind(dna.seq , pop=NULL, exp.char=c("a","t","g","c"), polyThres=1/100)

bin.dna@tab %>%
as.data.frame %>% 
.[,seq(from=1, to = dim(.)[2], by = 2 )] ->
simplified.bin.DNA

###
###
which(names(simplified.bin.DNA) =="1019.t")

#100:140

simplified.bin.DNA[-c(50:51),100:140] %>% 
PCA(graph = F) ->
pca.obj

pca.obj$ind$coord %>%
as.data.frame %>% 
mutate(ind = rownames(.)) %>%
#left_join(mutate(genotype, ind = paste(tip,".0", sep ="")  )) %>%  
ggplot(aes(
  x=Dim.2,
  y=Dim.3,
  #color = kar,
  #shape = pop
)) +
geom_point(size = 3) ->
pca.dna

ggsave(pca.dna, file = "pca.dna.pdf")

dimdesc(pca.obj) -> ddst
ddst$Dim.2 %>% as.data.frame %>% .[complete.cases(.),]

### pos 3844 in ALN

##
ggplot(tree, aes(x, y), 
branch.length="none"
) + 
geom_tree() + 
theme_tree() +
 geom_tiplab(size=2, align=TRUE, linesize=.5) + 
geom_nodelab(aes(x=branch, label=label), vjust=-.5, size=3) -> tree.plot


msaplot(tree.plot, dna.seq, 
window=c(998, 1008),  
offset = 5) -> dna.tree
ggsave(dna.tree, file = "dna.tree.tree.plot.pdf") 

q("no")

##
##
##



