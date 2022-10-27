##### ----> extract the In2Lt markers
##### 
##### 

### load modules
rm(list = ls())

library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(vroom)

### load
### 
matched.controls <- vroom("./Matched.controls.In2Lt.txt")

###
ag.folder <- "/project/berglandlab/DEST_fst_out/fst.ag.folder"

files.vec = system(paste("ls ",  ag.folder),
                   intern= T)

####
####
####
 
args = commandArgs(trailingOnly=TRUE)
i=as.numeric(args[1])

  message(paste(files.vec[i], sep = " | "))
  
  tmp.load = get(load(paste(ag.folder, files.vec[i], sep = "/")))
  
  selected.snps = c( unique(matched.controls$matched_to),
                     unique(matched.controls$control_snp)
                     )
  
  message("Now appliying the filter")
  tmp.load %>% 
    filter(SNP_id %in% selected.snps) -> 
    tmp.selected
  

  save(tmp.selected,
       file = paste("collect.job", i, "Rdata", sep = ".")
       )
