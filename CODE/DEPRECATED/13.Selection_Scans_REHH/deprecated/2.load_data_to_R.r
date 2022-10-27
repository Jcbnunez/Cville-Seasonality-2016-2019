### 
### Run EEH for CM samples
### 

library(rehh)
library(patchwork)
library(tidyverse)
library(magrittr)
library(ggbeeswarm)
library(data.table)
library(reshape2)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(viridis)
library(foreach)
library(doMC)
library(lubridate)
library(forcats)
library(viridis)
registerDoMC(2)

setwd("/scratch/yey2sn/Overwintering_ms/13.Selection_REHH")

vcf_file = "./VA_ch.HomoKarGenos.recode.vcf.gz"
hh <- data2haplohh(hap_file = vcf_file,
                   polarize_vcf = FALSE,
                   vcf_reader = "data.table")

save(hh, file = "rehh.object.CM.Rdata")