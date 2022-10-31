### creating enrichment plots
#use or and proportion from gwas/glm data to see relationship
#goal- make graph showing peaks of sliding window analysis accross chromsoomes and groups
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(cowplot)
#use doParallel package to register multiple cores that can be used to run loops in parallel
registerDoParallel(5)
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/May_2022_objects/")
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/March_2022_objects/")
#######################################################3
#gather together the files from the crescent plot analyss
###################################################

## import files
fl <- list.files("/scratch/bal7cg/Deficiency-Line-confirmation/gwas_glm_coenrich3/", pattern = "VA", full.names=T)
#only need first 101
#fl = fl[1:101]
groups = read_csv("/scratch/bal7cg/Deficiency-Line-confirmation/phenogroups3.22.csv")

groups = groups[,c(1,5)]
gwas.o <- foreach(fl.i=fl)%dopar%{
  #fl.i <- fl[2]
  message(fl.i)
  x=load(fl.i)
  gwas.o = peak.o
  #remake the inversion column to be descriptive
  gwas.o$inv.st = ifelse(gwas.o$inv == T, "inverted","non-inverted")
  #add fixed phenotype column
  gwas.o$phenotype = tstrsplit(x = gwas.o$gwas.pheno, split = ".", fixed = T)[[1]]
  #remove u's
  #gwas.o$phenotype = enc2utf8(gwas.o$phenotype)
  gwas.o$phenotype = gsub("μ","",gwas.o$phenotype)
  mergedata = merge(gwas.o, groups, by = "phenotype")
  mergedata$perm.st = ifelse(mergedata$glm.perm == 0, "observed","permutation")
  #we can use t.test to get confidence intervals for or and prop
  #but we need to do so while also dividing for threshold/chr/inv
  #ok- game plan is make a dataframe that has the ci, then merge it back onto mergedata
  
  #ok - now set chr/thr/inv as keys and merge the confidence intervals into there
  
  
  
  
  
  
  # gwas.win.o[,pa:=p.adjust(gwas.win.o$fet.p, method = "fdr")]
  
  mergedata
}
gwas.o <- rbindlist(gwas.o)
saveRDS(gwas.o, "/scratch/bal7cg/Deficiency-Line-confirmation/multi-peak.coenrich")

#load in data (from no grm gwas)

data = readRDS("multi-peak.coenrich")


#we only care about model F for now
gwas.o = data[glm.mod == 4]
#something went wrong with how the 95 % conficence intervals were calculated. we'll stick witht he old method for now, and delete those rows
gwas.o = gwas.o[,-c(4:7)]
#now- we want to estimate confidence intervals for each chr/thre/inv combination
# example = gwas.o[glm.perm == 1][chr == "2L"][inv == T]
# quantile(gwas.o$or, 0.05)
# can use group_by plus quantile
testdata = gwas.o %>%
  group_by(chr,inv, glm.perm, thr) %>%
  summarise( lowerlimitor = quantile(or, 0.01, na.rm = T),
             upperlimitor = quantile(or, 0.99, na.rm = T),
             lowerlimitprop = quantile(prop, 0.01, na.rm = T),
             upperlimitprop = quantile(prop, 0.99, na.rm = T)) %>%
  as.data.table(.)
#ok- now make a table of mean bounds, organized by thr, chr, and inv (only for permutations)
bounds = testdata %>% 
  filter(glm.perm != 0) %>% 
  group_by(chr, inv, thr) %>% 
  summarise(mean.lower.or = mean(lowerlimitor),
            mean.upper.or = mean(upperlimitor),
            mean.lower.prop = mean(lowerlimitprop),
            mean.upper.prop = mean(upperlimitprop)) %>% 
  as.data.table(.)

#for now we'll stick with the more concertavite threshold
bounds = bounds[thr == 0.05]
#confidence intervals- are they about the same for all permutations?
permdata = gwas.o[perm.st == "permutation"]
observedata = gwas.o[perm.st == "observed"][thr == 0.05]
save(bounds, observedata, file = "co.enrich.files" )
