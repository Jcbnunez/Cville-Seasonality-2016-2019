module load exonerate
module load samtools
module load bcftools

## position list

#### ---> look for motiff *TCCTGGAGATC* -- pos in exon 3835
#### ---> look for motiff *TCCTG-G-AGATC* -- pos in exon 3835

#### Official position list:
# --> Yakuba: NC_052527.2:9456177
# --> Sechelia: NC_045949.1:5099990
###step 1 find exonic region of yak to yak
query=Dyak.msp300.eq.txt
REFERENCE=/project/berglandlab/Dmel_genomic_resources/Other_drosophilids/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_genomic.fna
exonerate --model  est2genome  --maxintron 0 --percent 95 $query $REFERENCE > search.yakTOyak.txt

#         Query: D.yak.msp300.equ
#        Target: NC_052527.2 Drosophila yakuba strain Tai18E2 chromosome 2L, Prin_Dyak_Tai18E2_2.1, whole genome shotgun sequence
#         Model: est2genome
#     Raw score: 10005
#   Query range: 0 -> 2001
#  Target range: 9455179 -> 9457180

###step 2 find homology between regions in Mel abd Yak region of yak to yak
#samtools faidx $REFERENCE
samtools faidx $REFERENCE 'NC_052527.2:9455179-9457180' > Dyak.MSP300.ref.fasta

REFERENCE=Dmel.msp300.exact.txt 
query=Dyak.MSP300.ref.fasta
exonerate --model  est2genome  $query $REFERENCE > search.yakTomel.txt
### position in Yakuba = 9455179+999-1 = 9456177
#BAM=/scratch/yey2sn/Overwintering_ms/msp300.case/mapping_output/yak/SRR5860601/1_19.srt.rmdp.bam
#BAM=/scratch/yey2sn/Overwintering_ms/msp300.case/mapping_output/yak/SRR5860595/BIOKO_NE_4_6.srt.rmdp.bam
BAM=/scratch/yey2sn/Overwintering_ms/msp300.case/mapping_output/yak/SRR5860655/3_23.srt.rmdp.bam
samtools index $BAM

REFERENCE=/project/berglandlab/Dmel_genomic_resources/Other_drosophilids/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_genomic.fna
samtools tview $BAM $REFERENCE -p 'NC_052527.2:9456177'

####################################
####################################
####################################
####################################
### Sechelia ###
####################################
####################################
####################################
####################################


###step 1 find exonic region of sech to sech
query=Dsech.msp300.eq.txt
REFERENCE=/project/berglandlab/Dmel_genomic_resources/Other_drosophilids/GCF_004382195.2_Dsechelia_ASM438219v2_genomic.fna
exonerate --model  est2genome  --maxintron 0 --percent 95 $query $REFERENCE > search.sechTOsech.txt

#         Query: D.sech.msp300.equ
#        Target: NC_045949.1 Drosophila sechellia strain sech25 chromosome 2L, ASM438219v1, whole genome shotgun sequence
#         Model: est2genome
#     Raw score: 9990
#   Query range: 0 -> 1998
#  Target range: 5098992 -> 5100990

samtools faidx $REFERENCE 'NC_045949.1:5098992-5100990' > Dsech.MSP300.ref.fasta

REFERENCE=Dmel.msp300.exact.txt 
query=Dsech.MSP300.ref.fasta
exonerate --model  est2genome  $query $REFERENCE > search.SechTomel.txt
### position in Sechelia = 5098992+999-1 = 5099990

#BAM=/scratch/yey2sn/Overwintering_ms/msp300.case/mapping_output/sech/SRR5860633/DenisNF13.srt.rmdp.bam
BAM=/scratch/yey2sn/Overwintering_ms/msp300.case/mapping_output/sech/SRR5860626/DenisNoni10.srt.rmdp.bam
samtools index $BAM
REFERENCE=/project/berglandlab/Dmel_genomic_resources/Other_drosophilids/GCF_004382195.2_Dsechelia_ASM438219v2_genomic.fna
samtools tview $BAM $REFERENCE -p 'NC_045949.1:5099990'


####################################
####################################
####################################
####################################
### Mauritania ###
####################################
####################################
####################################
####################################

###step 1 find exonic region of sech to sech
query=Dsech.msp300.eq.txt
REFERENCE=/project/berglandlab/Dmel_genomic_resources/Other_drosophilids/GCF_004382145.1_Dmauritiana_ASM438214v1_genomic.fna
exonerate --model  est2genome  --maxintron 0 --percent 95 $query $REFERENCE > search.sechTOMau.txt

#         Query: D.sech.msp300.equ
#        Target: NC_046667.1 Drosophila mauritiana strain mau12 chromosome 2L, ASM438214v1, whole genome shotgun sequence
#         Model: est2genome
#     Raw score: 9495
#   Query range: 0 -> 1998
#  Target range: 5096167 -> 5098165

samtools faidx $REFERENCE 'NC_046667.1:5096167-5098165' > Dmau.MSP300.ref.fasta

REFERENCE=Dmel.msp300.exact.txt 
query=Dmau.MSP300.ref.fasta
exonerate --model  est2genome  $query $REFERENCE > search.MauTomel.txt
### position in Sechelia = 5096167+999-1 = 5097165
#samtools index $BAM

REFERENCE=/project/berglandlab/Dmel_genomic_resources/Other_drosophilids/GCF_004382145.1_Dmauritiana_ASM438214v1_genomic.fna
BAM=/scratch/yey2sn/Overwintering_ms/msp300.case/mapping_output/mau/SRR483621/*.RG.bam
samtools tview $BAM $REFERENCE -p 'NC_046667.1:5097165'



