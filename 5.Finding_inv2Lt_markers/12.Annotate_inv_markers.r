# Annotate the inversion markers.
# R script
# 
#module load intel/18.0 intelmpi/18.0
#module load goolf/7.1.0_3.1.4
#module load gdal proj R/4.0.0
#R

library(data.table)
library(tidyverse)
library(vcfR)

DGRP_inversion_status <- fread("/scratch/yey2sn/Overwintering_ms/Inversion_markers/inv_2lt_dgrp_lines.txt")

inversion_markers <- fread("/scratch/yey2sn/Overwintering_ms/Inversion_markers/in2lt_ld_47snps_informative_markers.txt", head = F)
names(inversion_markers) = "InvMarkers"

dgrp_vcf <- read.vcfR("/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz")

dgrp_vcf@fix %>% 
  as.data.frame ->
  dgrp_allelic_states
  
# Slice out markers
dgrp_allelic_states %>%
  filter(ID %in% inversion_markers$InvMarkers) ->
  inv_markers_states

#Slice out individuals
dgrp_vcf@gt %>% 
  as.data.frame ->
  dgrp_gts

#create gt filter
inv_markers_idRows <- which(dgrp_allelic_states$ID %in% inversion_markers$InvMarkers)
dgrp_gts[inv_markers_idRows,] ->
  GDRP_inverstion_gts

names(GDRP_inverstion_gts) = gsub("line", "DGRP", names(GDRP_inverstion_gts) )

DGRP_inversion_status %>%
  filter(`In(2L)t` %in% c("ST","INV")) %>%
  mutate(line_stat = paste(`DGRP Line`,`In(2L)t`, sep = "_")) -> inv_st_lines

GDRP_inverstion_gts[,inv_st_lines$`DGRP Line`] -> reduced_gt
names(reduced_gt) = inv_st_lines$line_stat

### Merge datasets
cbind(inv_markers_states, reduced_gt) -> inv_markers_states_gt

write.table(inv_markers_states_gt,
            file = "DGRP_inversion_states.txt", 
            append = FALSE, 
            quote = FALSE, 
            sep = " ",
            eol = "\n", 
            na = "NA", dec = ".", 
            row.names = FALSE,
            col.names = TRUE, 
            qmethod = c("escape", "double"),
            fileEncoding = "")

inv_markers_states_gt[,1:5] ->
  inversion_state_data

names(inversion_state_data)[4:5] = c("REF_ST", "ALT_INV")

write.table(inversion_state_data,
            file = "Inv2Lt_allele_states.txt", 
            append = FALSE, 
            quote = FALSE, 
            sep = " ",
            eol = "\n", 
            na = "NA", dec = ".", 
            row.names = FALSE,
            col.names = TRUE, 
            qmethod = c("escape", "double"),
            fileEncoding = "")




