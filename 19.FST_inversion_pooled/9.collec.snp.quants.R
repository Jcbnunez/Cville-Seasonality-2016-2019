### collect and calculate the quantile of SNPs
### 
### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(foreach)
library(doParallel)

snp.out.folder = "/scratch/yey2sn/Overwintering_ms/19.inv.fst/ind.snps.out"

files.vec = system(paste("ls ",  snp.out.folder),
                   intern= T)

args = commandArgs(trailingOnly=TRUE)
i=as.numeric(args[1])

                          tmp.snp <- get(load(  paste(snp.out.folder, files.vec[i], sep = "/")  ))
                          
                          anchor = str_split(files.vec[i], "\\.")[[1]][1]
                          
                          tmp.snp$year_diff[is.na(tmp.snp$year_diff)] = "space"
                          
                          tmp.snp$comp.set %>% unique() -> sets
                          tmp.snp$year_diff %>% unique() -> year_diffs.set
                          tmp.snp$pop1 %>% unique() -> pops.set
                          
                          #######
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
                          #######
                          
                          snp.metric %>%
                            .[complete.cases(.$perc.above.control),] ->
                            snp.metric.out
                            
                          save.folder = "/scratch/yey2sn/Overwintering_ms/19.inv.fst/quantile.ind.snps.out/"
                          
                          file.name = paste(save.folder, anchor, ".quantile.Rdata" , sep = "" )
                          
                          save(snp.metric.out, file = file.name)

