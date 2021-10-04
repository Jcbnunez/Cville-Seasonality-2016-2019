## Check file name pairings
library(tidyverse)
library(magrittr)

setwd("~/Dropbox/2021_Cville_DEST/Alyssa\'s\ Single\ Ind/samps_Mapping_pairs/")

all_mapping_pairs <- read.delim("~/Dropbox/2021_Cville_DEST/Alyssa\'s\ Single\ Ind/samps_Mapping_pairs/all_mapping_pairs.txt", header=FALSE)
names(all_mapping_pairs) = c("File1","File2")


GoodSamps_list <- read.table("~/Dropbox/2021_Cville_DEST/Alyssa\'s\ Single\ Ind/samps_Mapping_pairs/GoodSamps_list.txt", quote="\"", comment.char="")
names(GoodSamps_list) = "goodsamps"


#######

all_mapping_pairs %<>%
  separate(File1,
           remove = F,
           into = c("lane_1","sample_1"),
           sep = "_1_") %>% 
  separate(File2,
           remove = F,
           into = c("lane_2","sample_2"),
           sep = "_2_") %>%
  mutate(Merged_name = gsub(".fastq.gz", "", paste(lane_1, sample_1, sep = "_"))) %>%
  separate(sample_2, remove = F, into = c("sample_name","SL_extension"), sep = "\\_") 


## Apply good samps filter
all_mapping_pairs = 
  all_mapping_pairs %>%
  .[which(.$sample_name %in% GoodSamps_list$goodsamps),]

#Generate logical file
data.frame(
table(all_mapping_pairs$sample_1),
table(all_mapping_pairs$sample_2)
) %>%
  mutate(LOGIC_NAME = ifelse(Var1 == Var1.1, "yes","no"),
         LOGIC_NUMBER = ifelse(Freq == Freq.1, "yes8","no8")
         ) %>%
  mutate(unique_name = gsub(".fastq.gz", "" ,.$Var1) )-> logic_comparison_inds


### Save files
write.table(all_mapping_pairs,
            file = "Alyssa_ind_reads_guideFile.txt",
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t"
)

write.table(logic_comparison_inds,
            file = "Alyssa_ind_logic_comparison_inds.txt",
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t"
            )
