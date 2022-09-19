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
matched.controls <- vroom("./Matched.controls.quantiles.r075.In2Lt.txt")

###
matched.controls$matched_to %>% unique() -> in2lt_snps

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

save.folder = "/scratch/yey2sn/Overwintering_ms/19.inv.fst/ind.snps.out/"

file.name = paste(save.folder, snp.focus, ".snpdat.Rdata" , sep = "" )

save(tmp.selected, file = file.name)

