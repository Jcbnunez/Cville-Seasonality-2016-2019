##### ---> Correlated allele frequencies, pt 2, prepping for the matrix
##### 

library(tidyverse)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)

#####

load(file = "data_for_covariate_analysis_2L.Rdata")

# load in 
#dat_glm_p1, glm_2L_outP1perc, samps_metadata
#

## Selecting the sample pairs

sample_clusters <- c(
  "VA_ch_16_m07d08", #State 1
  "VA_ch_17_m07d07", #State 1
  "VA_ch_18_m07d12", #State 1
  
  "VA_ch_16_m09d02", #State 2
  "VA_ch_17_m08d17", #State 2
  "VA_ch_18_m09d06", #State 3
  
  "VA_ch_16_m09d16", #State 4
  "VA_ch_17_m10d12", #State 4
  "VA_ch_18_m09d20", #State 4
  
  "VA_ch_16_m10d14", #State 5
  "VA_ch_17_m10d26", #State 5
  "VA_ch_18_m10d18", #State 5
  
  "VA_ch_16_m11d11", #State 6
  "VA_ch_17_m11d09",  #State 6
  "VA_ch_18_m11d01" #State 6
)

### updating metadata
samps_metadata %>%
  filter(sampleId %in% sample_clusters) %>% 
  mutate(sample_cluster_cov = case_when(
    sampleId %in% c("VA_ch_16_m07d08","VA_ch_17_m07d07","VA_ch_18_m07d12") ~  "State_0",
    sampleId %in% c("VA_ch_16_m09d02", "VA_ch_17_m08d17","VA_ch_18_m09d06") ~  "State_1",
    sampleId %in% c("VA_ch_16_m09d16","VA_ch_17_m10d12","VA_ch_18_m09d20") ~  "State_2",
    sampleId %in% c("VA_ch_16_m10d14","VA_ch_17_m10d26","VA_ch_18_m10d18") ~  "State_3",
    sampleId %in% c("VA_ch_16_m11d11", "VA_ch_17_m11d09","VA_ch_18_m11d01") ~  "State_4"
  )) -> samps_metadata_working


samps_metadata_working$REP_id[grep("16",samps_metadata_working$sampleId)] = "REP_1"
samps_metadata_working$REP_id[grep("17",samps_metadata_working$sampleId)] = "REP_2"
samps_metadata_working$REP_id[grep("18",samps_metadata_working$sampleId)] = "REP_3"

####

glm_2L_outP1perc %>%
  filter(rnp.clean <= 0.01,
         invName == "2Lt") ->
  glm_2L_outP01perc

### generate combinations
unique_combs <- 
  combn(glm_2L_outP01perc$SNP_id ,2) %>% 
    t() %>%
  as.data.frame()

unique_combs %>%
  separate(V1, into = c("chr1", "pos1", "feature1"), remove = F ) %>%
  separate(V2, into = c("chr2", "pos2", "feature2"), remove = F ) %>%
  mutate(bp_dist = abs(as.numeric(pos1)-as.numeric(pos2))) ->
  unique_combs_sep
  
save(unique_combs_sep,
     glm_2L_outP01perc,
     dat_glm_p1,
     file = "data_for_step_3_covariance.Rdata")

