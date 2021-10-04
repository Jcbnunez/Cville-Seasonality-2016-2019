#!/usr/bin/env Rscript
# HEADER --------------------------------------------
#
# Author: JCB Nunez, PhD
# Copyright (c) JCB Nunez 2021
# Email:  yey2sn@virginia.edu
# 
# Date: May 14, 2021
#
# Script Name: Import GDS to R
#
# Script Description: This script imports DEST GDS into R
#
#
# Notes: Used the DEST 2.0 data
#
#

# allow user to pass inputs

args = commandArgs(trailingOnly=TRUE)

#Examples:

#Path to the GDS file
ingds=args[1]
#ingds="/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/PoolSNP_noRep_filter/dest.all.PoolSNP.001.50.10Mar2021.ann.noRep.gds"

#Path to the metadata
inmeta=args[2]
#inmeta="/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv"

#Name of the outfile to be appended here
outfile=args[3]
#outfile="DEST.2.0.poolSNP.Spatial.Temporal"

selected_chrs=c("2L", "2R", "3L", "3R")

# Load R packages

library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(tidyverse)
library(gmodels)
library(reshape2)
library(magrittr)

### open GDS file
  genofile <- seqOpen(ingds)

### get target populations
  samps <- fread(inmeta)

### Include DEST sets
  samps <- rbind(samps[set=="DrosRTEC"],
                samps[set=="DrosEU"],
                samps[set=="CvilleSet"]
                )

 ### get subsample of data to work on
  seqResetFilter(genofile)
  seqSetFilter(genofile, sample.id=samps$sampleId)

  snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        variant.id=seqGetData(genofile, "variant.id"),
                        nAlleles=seqNumAllele(genofile),
                        missing=seqMissing(genofile, .progress=T))

## choose number of alleles
snps.dt <- snps.dt[nAlleles==2]

seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt$variant.id)

snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

### select sites
  seqSetFilter(genofile, sample.id=samps$sampleId,
              snps.dt[chr%in%selected_chrs]$variant.id)
  
### get allele frequency data
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")

  dat <- ad$data/dp #if poolSNP
  #dat <- ad/dp #if SNAPE
  dim(dat)  

## Add metadata
colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")

#Add metadata ad
colnames(ad$data) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(ad$data) <- seqGetData(genofile, "sample.id")

#Add metadata dp
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")

#Generate metadata
left_join(data.frame(sampleId=rownames(dat)), as.data.frame(samps)) -> DEST_DGN_metadata


## Calculate effective coverage
### Generate a function to estimate effective coverage
### Nc = (1/N + 1/R) - 1
calc_effcov = function(NomCOV, FlyPool) {
  
  EC1 = ((2*FlyPool)*NomCOV)
  EC2 = ((2*FlyPool)+NomCOV)
  Nc=(EC1/EC2)
  return(Nc)
  
}

#######
left_join(data.frame(sampleId=rownames(dp)), as.data.frame(samps)) -> DP_DEST_DGN_metadata

dp %>%
  .[,sample(dim(.)[2], 100000 )] %>%
  as.data.table() %>%
  mutate(sampleId = rownames(dp)) %>%
  left_join(DP_DEST_DGN_metadata[,c("sampleId","nFlies")]) %>%
  melt(id = c("sampleId","nFlies")) %>% 
  mutate(EFFCOV = calc_effcov(value, nFlies) ) %>% 
  group_by(sampleId) %>%
  summarise(MeanEC = mean(EFFCOV, na.rm = T)) -> EFFCOV_samps

#Study of effective coverage
left_join(samps, EFFCOV_samps ) -> samps_EFFCOV

samps_EFFCOV$city = gsub("Charlotttesville","Charlottesville", samps_EFFCOV$city)
samps_EFFCOV$city = gsub("Odesa","Odessa", samps_EFFCOV$city )
samps_EFFCOV$city[grep("Yesiloz", samps_EFFCOV$city )] = "Yesiloz"
samps_EFFCOV$city[grep("Chornobyl", samps_EFFCOV$city )] = "Chernobyl"
samps_EFFCOV$city[grep("Kyiv", samps_EFFCOV$city )] = "Kyiv"

samps_EFFCOV %>%
  group_by(set) %>%
  summarise(N=n())

samps_EFFCOV %>%
  group_by(country) %>%
  summarise(N=n()) %>% as.data.frame()

samps_EFFCOV %>%
  group_by(year) %>%
  summarise(N=n()) %>% as.data.frame()


samps_EFFCOV$MeanEC %>% 
  .[which(samps_EFFCOV$set == "CvilleSet")] %>%    
  ci %>%
  round(3)

samps_EFFCOV$MeanEC %>% 
  .[which(samps_EFFCOV$set == "CvilleSet")] %>%    
  quantile() %>%
  round(3)


## Generate filter EFFCOV to > X= for charlottesville
  samps_EFFCOV %>%
    .[which(samps_EFFCOV$set == "CvilleSet")] %>%    
    .[which(.$MeanEC < 30),] %>%
    .$sampleId ->
    Cville_to_remove

  samps_EFFCOV %>%
    .[-which(.$sampleId %in% Cville_to_remove),] ->
    samps_filter_EC
  
  
  samps_filter_EC %>%
    group_by(set) %>%
    summarise(N = n()) 

  samps_filter_EC %>%
    group_by(city, year) %>%
    summarise(N = n()) %>%
    dcast(city~year) -> samps_year
  
  samps_year %<>%
    mutate(samps_per_year = rowMeans(samps_year[,-1], na.rm = T),
           NA_count= apply(., 1, function(x) sum(is.na(x))) )
  
  samps_year %>%
    .[which(.$samps_per_year >= 2  &
            .$NA_count < 10),] ->
    samps_year_to_use
 
## Make figure:

  samps_EFFCOV %>%
    .[which(.$city %in% samps_year_to_use$city),] %>%
    .[-which(.$year == 2012 & .$city == "Charlottesville"),]  -> 
    filtered_samps_for_analysis

  filtered_samps_for_analysis %>%
    .[which(.$city %in% samps_year_to_use$city),] %>%
    group_by(country) %>%
    summarise(N= n())
  
  filtered_samps_for_analysis %>%
    .[which(.$city %in% samps_year_to_use$city),] %>%
    group_by(country, city) %>%
    summarise(N= n())
  
  filtered_samps_for_analysis %<>% 
    mutate(Month = 
             month(as.Date(collectionDate, 
                                               format = "%m/%d/%Y")))
  
    
  filtered_samps_for_analysis %>% 
    ggplot(
      aes(
        y=paste(country, city, sep = " "),
        x=as.numeric(Month),
        fill = factor(month.abb[Month], levels = month.abb),
      )
    ) +
    ylab("Locale") +
    xlab("Month") +
    geom_point(size = 3.2,
               shape = 21) + 
    scale_fill_brewer(name = "Month",
                      palette = "Spectral") +
    xlim(5.5,12.5) +
    theme_bw() +
    facet_wrap(~year, nrow = 1)->
    Sample_schedule
  
  ggsave(Sample_schedule, 
         file = "Sample_schedule.pdf",
         width = 12,
         height = 4)
  ggsave(Sample_schedule, 
         file = "Sample_schedule.png",
         width = 12,
         height = 4)

## Apply filter EFFCOV to > X+
dat %>% 
  .[which(rownames(.) %in% filtered_samps_for_analysis$sampleId),] %>%
  t() %>%
  as.data.frame() %>%
  mutate(SNP_id = rownames(.)) -> dat_filtered_t

dat_filtered_t %>%
  dim

### Save object
save(filtered_samps_for_analysis,
	dat_filtered_t,
	file = paste(outfile,"ECfiltered.Rdata", sep = "."))
