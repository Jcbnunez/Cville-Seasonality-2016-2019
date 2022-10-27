### This script prepares data for haplovalidate -- pt 2

library(tidyverse)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)


load("./sync_format_table_interm_file.Rdata")
load("./ad_dp_obj.Rdata")

#####
#####

sync_format_table %>%
  dcast(chr+pos~variable, value.var = "sync_form") ->
  sync_format_table_dcast


clusters <- c(
  "VA_ch_16_m07d08", #State 1
  "VA_ch_17_m07d07", #State 1
  "VA_ch_18_m07d12", #State 1
  
  "VA_ch_16_m09d02", #State 2
  "VA_ch_17_m08d17", #State 2
  "VA_ch_18_m09d06", #State 3
  
  "VA_ch_16_m09d16", #State 4
  "VA_ch_17_m10d12", #State 4
  "VA_ch_18_m09d20", #State 4
  
  "VA_ch_16_m10d14", #State 5
  "VA_ch_17_m10d26", #State 5
  "VA_ch_18_m10d18", #State 5
  
  "VA_ch_16_m11d11", #State 6
  "VA_ch_17_m11d09",  #State 6
  "VA_ch_18_m11d01" #State 6
)

sync_format_table_dcast_df <- as.data.frame(sync_format_table_dcast)

left_join(ad_dp_obj[,c("chr", "pos", "REF")], sync_format_table_dcast_df[,c("chr",  "pos" ,clusters)] ) ->
  sync_file

sync_file %>% head

#save a text file
fwrite(sync_file, 
       file = "DEST_Cville_sync.selectedClusts.2L.sync", 
       append = FALSE, 
       quote = FALSE, 
       sep = "\t",
       row.names = FALSE,
       col.names = FALSE)

