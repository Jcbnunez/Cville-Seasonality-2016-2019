
### Run PCA on DEST data
### 
rm(list = ls())

#load data
# This R object was premade in script 1 ==> 1.Import_GDStoR.r
data_in="/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST.2.0.poolSNP.Spatial.Temporal.mAF_Miss_Mean_Filt.ECfiltered.Rdata"
name="all"

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

### load data
load(data_in)

### Run up to here to begin

dat_for_Analysis %>%
  as.data.frame() %>% 
  PCA(scale.unit = F, graph = F, ncp = 20) ->
  PCA_object

save(PCA_object, 
     file = "/scratch/yey2sn/Overwintering_ms/2.Temporal_Spatial_structure/PCA.object.all.Rdata")

