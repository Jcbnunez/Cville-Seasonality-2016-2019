### R Script for talk and paper
### 
## This script makes figure 3 of the paper

rm(list = ls())
# Load packages

#args = commandArgs(trailingOnly=TRUE)

## the 2 arguments
## Argument 1 is the Rdata with the genotype matrix
## Argument 2 is the city to analyze
## Argument 3 is the name of the analysis

#data_in=args[1]
data_in="/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST.2.0.poolSNP.Spatial.Temporal.AllDat.ECfiltered.Rdata"

#city_target=unlist(strsplit(args[2], ",") )

#name=args[3]
name="global_set"

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
library(scales)
library(gmodels)

##Import the inversion mapping SNPs
inversions = read.table(
  "/project/berglandlab/Dmel_genomic_resources/Inversions/inversion_makers_kapun.txt",
  head = T)

names(inversions)[2:3] = c("chr","pos")

## Load the genotype matrix
load(data_in)

#### Begin Measure inversions
### Measure inversion proportions
inversions %>% head

dat_filtered_t %>% 
  separate(SNP_id, 
           into = c("chr","pos"), 
           sep = "_", 
           remove = F) ->
  dat_filtered_t_chr_pos

dat_filtered_t_chr_pos$pos = as.numeric(dat_filtered_t_chr_pos$pos)

left_join(inversions, 
          dat_filtered_t_chr_pos, 
          by = c("chr","pos")) %>%
  .[complete.cases(.),] %>%
  melt(id = c("inversion", 
              "chr", 
              "pos", 
              "SNP_id", 
              "allele")) %>% 
  group_by(inversion,variable) %>%
  summarize(Inv_freq = mean(value)) %>% 
  dcast(variable~inversion, value.var = "Inv_freq") ->
  Inversion_frequencies

names(Inversion_frequencies)[1] = "sampleId" 

left_join(filtered_samps_for_analysis, Inversion_frequencies) ->
  filtered_samps_for_analysis
########## <--- End measure inverions

###
###
###

dat_filtered_t %>% ###<-- this object was created above
  .[,-which(names(.) %in% c("SNP_id"))] %>%
  t() ->
  dat_AF 

filtered_samps_for_analysis$city = gsub("Charlotttesville","Charlottesville", 
                                        filtered_samps_for_analysis$city )

filtered_samps_for_analysis$city = gsub("Odesa","Odessa", 
                                        filtered_samps_for_analysis$city )

##filtered_samps_for_analysis %>%
##  .[which(.$city %in% city_target ),] %>%
##  .$sampleId -> samps_target
##
##if("Charlottesville" %in% city_target) {
##  samps_target[-which(samps_target
##                      %in% c("VA_ch_12_spring", 
##                             "VA_ch_12_fall" ))] ->
##    samps_target
##}

dat_AF %>%
  t() %>% 
  as.data.frame -> dat_AF_samps_target

## Some characterizations of AFs and subsequent filtering
MeanAF=c()
MinAF=c()

apply(dat_AF_samps_target,
      1, FUN=mean, na.rm=TRUE ) -> MeanAF
data.frame(SNP_id = dat_filtered_t$SNP_id, MeanAF) -> MeanAF

apply(dat_AF_samps_target,
      1, FUN=min, na.rm=TRUE ) -> MinAF
data.frame(SNP_id = dat_filtered_t$SNP_id, MinAF) -> MinAF

cbind(dat_AF_samps_target, MeanAF, MinAF[-1]) -> dat_AF_samps_target

##
dat_AF_samps_target %>%
  .[which(.$MeanAF > 0.00 & .$MeanAF < 1.00),] %>%
  .[which(.$MinAF > 0.001),] ->  ### This samples only polymorphic sites
  dat_AF_samps_target_filtered

###
count_NA = function(x){
  return(sum(is.na(x)))
}
MissDat=c()

pool_cols = grep("SNP_id|MeanAF|MinAF", colnames(dat_AF_samps_target_filtered), invert = T)

apply(dat_AF_samps_target_filtered[,pool_cols],
      1, FUN=count_NA ) -> MissDat

n_pools = length(pool_cols)

data.frame(SNP_id = dat_AF_samps_target_filtered$SNP_id, missing_rate = c(MissDat/n_pools) ) -> MissDat

cbind(dat_AF_samps_target_filtered, MissDat[-1]) -> dat_AF_samps_target_filtered

dat_AF_samps_target_filtered %>%
  filter(missing_rate < 0.01) ->  ### This samples only polymorphic sites
  dat_AF_samps_target_filtered

###

dat_AF_samps_target_filtered %>%
  separate(SNP_id, into = c("chr","pos"), sep = "_") %>% 
  summarise(N=n())

dat_AF_samps_target_filtered %>%
  separate(SNP_id, into = c("chr","pos"), sep = "_") %>% 
  group_by(chr) %>%
  summarise(N=n())

# Remove X
dat_AF_samps_target_filtered %>%
  .[grep("X", .$SNP_id, invert = T ),] %>% 
  .[,-which(names(.) %in% c("SNP_id", "MeanAF", "MinAF", "missing_rate"))] %>% 
  t() ->
  dat_for_Analysis

colnames(dat_for_Analysis) = dat_AF_samps_target_filtered$SNP_id[grep("X", dat_AF_samps_target_filtered$SNP_id, invert = T )]

save(dat_for_Analysis, 
     file = "/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST.2.0.poolSNP.Spatial.Temporal.mAF_Miss_Mean_Filt.ECfiltered.Rdata")
