### make PCA from Cville chromosomes 
### 

rm(list = ls())
# Load packages

#load packages
library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(gtools)
library(patchwork)
library(ggalt)
library(reshape2)
library(ggrepel)
library(sp)
library(permute)
library(RColorBrewer)
library(devtools)
library(lubridate)


######

objects <- c(
  "/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_2L.ECfiltered.Rdata",
  "/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_2R.ECfiltered.Rdata",
  "/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_3L.ECfiltered.Rdata",
  "/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_3R.ECfiltered.Rdata"
)

pca_table_list = list()
for(i in 1:length(objects)){
  
  load(objects[i])
  o %>% colnames() %>% .[1] -> lead_snp
  
  lead_snp %>% data.frame(headsnp = .) %>% 
    separate(headsnp , into = c("CHR","POS")) -> lead_snp_guide
  
  print(lead_snp_guide$CHR)
  
  filtered_samps_for_analysis %>%
    filter(city == "Charlottesville",
           MeanEC > 30) %>% 
    .$sampleId -> select_samples
  
  o %>%
    as.data.frame() %>% 
    filter(rownames(.) %in%  select_samples) ->
    snps.tmp
    
  pca_table_list[[i]] = snps.tmp
  
}

chrs_o = do.call(cbind, pca_table_list)


