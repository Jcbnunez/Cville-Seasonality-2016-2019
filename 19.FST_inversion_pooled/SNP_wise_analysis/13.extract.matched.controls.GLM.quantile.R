##### ----> extract the In2Lt markers
##### 
##### 

### load modules
rm(list = ls())

library(data.table)
library(tidyverse)
library(magrittr)
library(foreach)
library(doParallel)
library(vroom)

args = commandArgs(trailingOnly=TRUE)
i=as.numeric(args[1])


### load
### 
matched.controls <- vroom("./Matched.controls.quantiles.GLM.txt")

###
matched.controls$matched_to %>% unique() -> in2lt_snps
### --> 3261

snp.focus = in2lt_snps[i]

###
ag.folder <- "/project/berglandlab/DEST_fst_out/fst.ag.folder"

files.vec = system(paste("ls ",  ag.folder),
                   intern= T)

####
####
####

tmp.selected =
foreach(k=1:length(files.vec),
        .combine = "rbind")%do%{
        
         message(paste(k ,files.vec[k], sep = " | "))
        
         tmp.load = get(load(paste(ag.folder, files.vec[k], sep = "/")))
         
         control.sets = matched.controls %>%
           filter(matched_to == snp.focus) %>%
           .$control_snp -> controls
         
         selected.snps = c( snp.focus, controls)
         
         message("Now appliying the filter")
         tmp.load %>% 
           filter(SNP_id %in% selected.snps) -> 
           tmp.selected.subset
         
         tmp.selected.subset %<>%
           mutate(anchor.set = case_when(SNP_id == snp.focus ~ "focus",
                                         SNP_id != snp.focus ~ "control"))
         
         return(tmp.selected.subset)
        }

#### --- save tmp subset to folder

tmp.snp <- tmp.selected

anchor = snp.focus

tmp.snp$year_diff[is.na(tmp.snp$year_diff)] = "space"

tmp.snp$comp.set %>% unique() -> sets
tmp.snp$year_diff %>% unique() -> year_diffs.set
tmp.snp$pop1 %>% unique() -> pops.set

snp.metric <- foreach(sets.i = sets, .combine = "rbind", .errorhandling = "remove"  )%do%{
  foreach(years.i = year_diffs.set, .combine = "rbind", .errorhandling = "remove" )%do%{
    foreach(pops.i = pops.set, .combine = "rbind", .errorhandling = "remove" )%do%{
      
      message(paste(sets.i, years.i, pops.i, sep= " | "))
      
      tmp.snp %>%
        filter(comp.set == sets.i & year_diff == years.i & pop1 == pops.i) %>%
        filter(anchor.set == "focus") ->
        tmp.snp.subset.anchor
      
      tmp.snp %>%
        filter(comp.set == sets.i & year_diff == years.i & pop1 == pops.i) %>%
        filter(anchor.set == "control") ->
        tmp.snp.subset.control
      
      sum(tmp.snp.subset.control$mean.SNPwise.FST < tmp.snp.subset.anchor$mean.SNPwise.FST) -> 
        percentile
      
      inner.out = data.frame(
        anchor.snp = anchor,
        comp.set = sets.i,
        year_diff = years.i,
        pop = pops.i,
        perc.above.control = percentile/length(tmp.snp.subset.control$mean.SNPwise.FST)
      )
      
      return(inner.out)
      
    }}}## foreach sets

####
snp.metric %>%
  .[complete.cases(.$perc.above.control),] ->
  snp.metric.out

save.folder = "/scratch/yey2sn/Overwintering_ms/19.inv.fst/quantile.ind.snps.out.glm/"

file.name = paste(save.folder, anchor, ".quantile.Rdata" , sep = "" )

save(snp.metric.out, file = file.name)

