### Map DGRP 2lt markers to Dmel6

library(data.table)
library(tidyverse)
library(magrittr)
library(reshape2)

#read in base data

liftOverFile <- "/project/berglandlab/Dmel_genomic_resources/liftOver_files/dest.all.PoolSNP.001.50.dm3.dm6.csv"
marker_2lt_file <- "/project/berglandlab/Dmel_genomic_resources/Inversions/DGRP_2lt_Markers/inv2L_informative_markers_Dm3.txt"


liftOver_marker <- fread(liftOverFile)
inversion_markers <- fread(marker_2lt_file)

#modify marker file
names(inversion_markers)[3:5] = c("SNP_id_dm3",
                                  "dm3_chr",
                                  "dm3_pos"  )

#merge files
#

left_join(inversion_markers, liftOver_marker) %>% 
  mutate(SNP_id_dm6 = paste(dm6_chr, dm6_pos, "SNP", sep = "_" )) ->
  updated_marker_list

write.table(updated_marker_list, 
            file = "/project/berglandlab/Dmel_genomic_resources/Inversions/DGRP_2lt_Markers/inv2L_informative_markers_Dm3toDm6.txt", 
            append = F, 
            quote = F, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = F,
            col.names = T)
