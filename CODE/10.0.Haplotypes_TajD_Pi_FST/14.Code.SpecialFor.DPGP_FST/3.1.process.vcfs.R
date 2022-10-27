### Clean up VCFs
### 

library(tidyverse)
library(vroom)

args = commandArgs(trailingOnly=TRUE)

file.o=as.character(args[1])

#file.o="/project/berglandlab/jcbnunez/emergency_scratch/dpgp3/vcfs_o/ZI104_Chr2L.seq.named.fasta.vcf"
  
vcf.r = vroom(file.o, delim = "\t", comment = "##")

samp.id = names(vcf.r)[10]
names(vcf.r)[10]= "samp.tmp"

### repair vcf
### 
vcf.r %>%
  mutate(samp.tmp = case_when(ALT == "N" ~ "./.",
                              REF != ALT ~ "1/1",
                              )) %>%
  mutate(ID = paste(POS,"SNP", sep = "_")) ->
  vcf.r.fix

vcf.r.fix %>%
  filter(ALT != "N") -> vcf.r.fix.r2

names(vcf.r.fix.r2)[10] = samp.id

vroom_write(vcf.r.fix.r2, paste("vcfs_diploidized/", samp.id,"R.fix",".vcf", sep = "") )

