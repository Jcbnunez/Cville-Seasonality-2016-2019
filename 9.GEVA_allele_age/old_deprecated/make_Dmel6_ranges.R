library(tidyverse)
library(magrittr)

#This script generates a range of chromosome intervals to run various programs on.
# by Joaquin Nunez

chr_names =c("X", "2L", "2R", "3L", "3R")
chr_Ls = c(23542271, 23513712, 25286936, 28110227, 32079331)
data.frame(chr_names, chr_Ls) -> Dmel6_metadat


## start i loop
allchrs_list = list()
for(i in 1:dim(Dmel6_metadat)[1]){

tmp_chr = Dmel6_metadat[i,"chr_names"]
tmp_ranges = seq(from = 1, 
    to = Dmel6_metadat[i,"chr_Ls"],
    by= 100000)

tmp_ranges[which(tmp_ranges == tail(tmp_ranges, n=1))] = Dmel6_metadat[i,"chr_Ls"]

tmp_size = length(tmp_ranges)
tmp_list = list()

  for(j in 1:tmp_size) {
    print(j)
    if(j == tmp_size){ print( "im done" )}
    else {
      out_df = data.frame(start = NA, end= NA, chr = NA)
      out_df$start[1] = tmp_ranges[j]
      out_df$end[1] = tmp_ranges[j+1]
      out_df$chr[1] = tmp_chr
      
      tmp_list[[j]] = out_df
      
      do.call(rbind,tmp_list) -> coordinates_tmp
    }#if
  }# j
  allchrs_list[[i]] = coordinates_tmp
}# i

  do.call(rbind,allchrs_list) -> guide_file_Dmel6

  guide_file_Dmel6 %<>%
    mutate(range = end - start,
           end_adj=format(as.numeric(end-1), scientific=F),
           start_adj=format(as.numeric(start), scientific=F),
           joint = paste(chr,":",
                                format(start_adj, scientific=F),
                                "-", 
                                format(end_adj, scientific=F), 
                                sep =""))
  
guide_file_Dmel6$joint = gsub("\ +", "", guide_file_Dmel6$joint)
guide_file_Dmel6$end_adj = gsub("\ +", "", guide_file_Dmel6$end_adj)
guide_file_Dmel6$start_adj = gsub("\ +", "", guide_file_Dmel6$start_adj)

write.table(guide_file_Dmel6,
            file = "array_job_guidefile.txt",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = F)
            
  
  
