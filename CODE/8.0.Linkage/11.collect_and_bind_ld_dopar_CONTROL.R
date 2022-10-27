### Code to process ld data
### 

library(tidyverse)
library(data.table)
library(magrittr)
library(forcats)
#parallel computing in R
library(foreach)
library(doMC)
registerDoMC(10) ## using an i-job
###
###
system(intern = T, "ls ld_control_plink | grep -v 'log' | grep -v 'nosex'") -> ld_files
#ld_list = list()

ld_list <- foreach(i=1:length(ld_files), 
                   #.combine = "rbind",
                   .errorhandling = "remove")%dopar%{
                     
                     print(i)
                     
                     tmp <- fread(paste("./ld_control_plink/", ld_files[i], sep = ""))
                     
                     mins <- apply(tmp[,c("BP_A","BP_B")], 1, min)
                     maxs <- apply(tmp[,c("BP_A","BP_B")], 1, max)
                     
                     tmp %<>%
                       mutate(pair_id = paste(mins,maxs, sep = "_" ) )
                     
                     return(tmp)
                     #ld_list[[i]] = tmp
                   }

ld_df = do.call(rbind, ld_list)

save(ld_df, file = "MATCHEDCONTROLs.merged.ld.Rdata")

