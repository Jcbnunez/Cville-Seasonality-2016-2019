## plot LD correlation analysis
## 

#Load libraries
library(tidyverse)
library(vcfR)
library(data.table)
library(adegenet)
library(FactoMineR)
library(magrittr)
library(reshape2)
library(zoo)


#load in the ld data

ld_in_files = list()
read_in_folder <- "/project/berglandlab/Dmel_genomic_resources/Inversions/DGRP_2lt_Markers/Rank90_99"

files <- system( paste("ls", read_in_folder, sep = " " ) , intern = T)

for(i in 1:length(files)){
  
  print(i/length(files)*100)
  
  tmp <- fread(paste(read_in_folder, files[i], sep = "/"), sep = "\ " )
  names(tmp)[1] = c("id_tag_rank")
  
  tmp$id_tag_rank =  gsub("\\t", "", tmp$id_tag_rank)
    
  tmp %<>%
    separate(id_tag_rank, into = c("id_tag", "rank"), sep = "\\|")
  
  ld_in_files[[i]] = tmp
  
}

ld_df = do.call(rbind, ld_in_files )

save(ld_df, file = "ld_df.rank90_99.Rdata")


