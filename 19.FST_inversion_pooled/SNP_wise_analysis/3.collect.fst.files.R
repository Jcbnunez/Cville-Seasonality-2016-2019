#### Collect files from the FST analysis
#### 
#### 

### load modules
rm(list = ls())

library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)

### data input
args = commandArgs(trailingOnly=TRUE)
chr.sel=as.character(args[1])
#chr.sel="2L"
  
root.folder <- "/project/berglandlab/DEST_fst_out/out.files"

files.vec = system(paste("ls ",  root.folder, "| grep 'SNPFST' ", " | grep", chr.sel),
                   intern= T)

merged.dfs = foreach(i = 1:length(files.vec),
                     .combine = "rbind")%do%{
  file.tmp = files.vec[i]
  message(paste(i, file.tmp, sep = " | "))
  
  tmp.load = get(load(paste(root.folder, file.tmp, sep = "/")))
  
  tmp.load %>%
    filter(comp.set %in% c( "time.cville","time.eu.e","time.eu.w", "core20.seasonality" ) ) %>%
    group_by(comp.set, SNP_id, year_diff, chr, pos, invName, same.locale, pop1) %>%
    summarise(
      mean.SNPwise.FST = mean(SNPwise.FST, na.rm = T),
      sd.SNPwise.FST = sd(SNPwise.FST, na.rm = T),
      median.SNPwise.FST = median(SNPwise.FST, na.rm = T),
    )  %>% mutate(analysis.type = "time") -> tmp.load.time.ag
  
  tmp.load %>%
    filter(comp.set %in% c( "space.eu.e","space.eu.w", "space.NoA.E" ) ) %>%
    group_by(comp.set, SNP_id, chr, pos, invName, same.locale) %>%
    summarise(
      mean.SNPwise.FST = mean(SNPwise.FST, na.rm = T),
      sd.SNPwise.FST = sd(SNPwise.FST, na.rm = T),
      median.SNPwise.FST = median(SNPwise.FST, na.rm = T),
    ) %>% mutate(year_diff = NA, analysis.type = "space", pop1 = comp.set)  -> tmp.load.space.ag
  
  tmp.load.ag.sp.tim = rbind(tmp.load.time.ag, tmp.load.space.ag)
  
  return(tmp.load.ag.sp.tim)
}


### Create sets to save

set.vec = c("core20.seasonality","space.eu.e","space.eu.w","space.NoA.E",
            "time.cville","time.eu.e","time.eu.w")

foreach(k = 1:length(set.vec))%do%{
  
  merged.dfs %>%
    filter(comp.set == set.vec[k]) ->
    flt.df
  
  save(flt.df,
       file = paste("SNPFST.AG", chr.sel, set.vec[k], "Rdata", sep = "." ) )

}
