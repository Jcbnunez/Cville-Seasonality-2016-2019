## haplovalidate
library(haploReconstruct)  
library(psych)
library(stringr)
library(data.table)
library(tidyverse)
library(haplovalidate)
###

## column names should be chr, pos and score 
load("/project/berglandlab/summarized_dest_glm/glm.out.VA_ch_0.Rdata")

glm.out %>% 
  filter(
         #rnp.clean<0.01,
         mod == "aveTemp",
         chr == "2L") ->
  glm_outP1perc

glm_outP1perc %>%
  select(chr, pos, b) -> out_snps

names(out_snps)[3] = c("score")

out_snps <- out_snps[order(out_snps$pos),]
out_snps$score = as.numeric(out_snps$score)

### load sync
### Extract outlier loci from sync file
sync_file <- fread("./DEST_Cville_sync.selectedClusts.2L.sync")

sync_file %>% 
  filter(V1 == "2L") -> 
  sync_file_2l

sync_file_2l %>% 
  group_by(V1, V2) %>% 
  slice_head() ->
  sync_2L_dedupl

#sync_2L_dedupl %>%
#  filter(V2 %in% out_snps$pos)  -> 
#  sync_glm
#sync_glm %>% dim

fwrite(sync_file_2l, 
       file = "DEST_Cville_sync.selectedClusts.2L.dedup.sync", 
       append = FALSE, 
       quote = FALSE, 
       sep = "\t",
       row.names = FALSE,
       col.names = FALSE)

######
syncfile_g <- "./DEST_Cville_sync.selectedClusts.2L.dedup.sync"

repl <- 1:3
gens <- c(0,3,6,9,12)

### define which columns of the sync file contain the base population
base.pops <- c(rep(TRUE, length(repl)),rep(FALSE,length(repl)*(length(gens)-1)))

### define which columns should be used to polarize for the rising allele (e.g. list(c(F0Rep1,F59Rep1),c(F0Rep2,F59Rep2),...))
polaRise = list(c(1,13),c(2,14),c(3,15)) 

### load frequency data
cands.all <- sync_to_frequencies(syncfile_g,
                                 base.pops=base.pops,
                                 header=FALSE,
                                 mincov=15,
                                 polaRise = polaRise)

cands.all %>% 
  .[complete.cases(.),] ->
  comp_cases

comp_cases %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) ->
  comp_cases_id

comp_cases_id %>% 
  group_by(SNP_id) %>% 
  slice_head() ->
  comp_cases_id_dedupl

comp_cases_id_dedupl %>%
  .[order(.$pos),] ->
  comp_cases_id_dedupl_sort

save(comp_cases_id_dedupl_sort, file = "cands.all.rdata")


