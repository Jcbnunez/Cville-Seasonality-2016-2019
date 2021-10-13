### finding SNPs in DGRP 2LT present on CM
### also do a sanity check
### 

# LOAD BCFTOOLS aourside of R
#load libraries
library(data.table)
library(tidyverse)
library(magrittr)

#import data
markers_metadata <- "/project/berglandlab/Dmel_genomic_resources/Inversions/DGRP_2lt_Markers/inv2L_informative_markers_Dm3toDm6.txt"
DGRP_vcf <- "/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.vcf.gz"
CM_vcf <- "/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.wSNPids.vcf.gz"
  
#make vcf data
system( paste("bcftools query  --regions 2L -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n'",
              DGRP_vcf,
              "> DGRP_vcf.snpids.2L.txt",
              sep = " "))

system( paste("bcftools query  --regions 2L -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n'",
              CM_vcf,
              "> CM_vcf.snpids.2L.txt",
              sep = " "))

#Load files
In2LT_markers = fread(markers_metadata)

DGRP_SNPs = fread("./DGRP_vcf.snpids.2L.txt")
names(DGRP_SNPs) = c("dm3_chr",  "dm3_pos", "dm3_REF", "dm3_ALT", "dm3_ID")

CM_SNPs = fread("./CM_vcf.snpids.2L.txt")
names(CM_SNPs) = c("dm6_chr",  "dm6_pos", "dm6_REF", "dm6_ALT", "dm6_ID")

#start comparison
In2LT_markers %<>%
  mutate(logical_lifOver = SNP_id_dm3 ==  SNP_id_dm6) 

In2LT_markers %>% 
  filter(logical_lifOver != FALSE) ->
  In2LT_markers_flt


left_join(In2LT_markers_flt, DGRP_SNPs) %>%
  left_join(CM_SNPs) %>% 
  filter(dm3_ID == dm6_ID) %>% 
  filter(dm3_REF == dm6_REF)  %>% 
  filter(dm3_ALT == dm6_ALT) ->
  In2LT_markers_flt_inSilicoValid

write.table(In2LT_markers_flt_inSilicoValid, 
            file = "Inv2L_markers_to_use_CM.txt", 
            append = F, 
            quote = F, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = F,
            col.names = T)

