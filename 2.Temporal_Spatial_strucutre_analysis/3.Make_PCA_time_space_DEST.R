
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
library(data.table)
library(foreach)
library(scales)

### load data
load(data_in)

### Prepare chromsome samples
extra_data <- "/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_2L.ECfiltered.Rdata"
load(extra_data)
head(filtered_samps_for_analysis)
###
## PLOT NEC
load("/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST_EC_metadata.Rdata")

samps_EFFCOV %>%
  filter(locality == "VA_ch") %>%
  ggplot(aes(
    x=sampleId,
    y=MeanEC,
    label=round(MeanEC,2),
    fill = year
  )) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, nudge_y = 5) +
  geom_hline(yintercept = 28) +
  ggtitle("Eff Cov tresh = 28X") +
  coord_flip() ->
  EC_barplot

ggsave(EC_barplot , file ="EC_barplot.pdf")

###
###
load("/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST_EC_metadata.Rdata")
samps_EFFCOV %<>%
  #filter(!(locality == "PA_li" & year <= 2012)) %>%
  filter(MeanEC > 30) ### apply filter here

samps_EFFCOV %>%
  filter(locality %in%
           c("DE_Bro","DE_Mun","FI_Aka","PA_li","TR_Yes","UA_Ode", "UA_od","VA_ch","WI_cp")) %>%
  group_by(locality, year) %>%
  summarize(N=n()) %>%
  dcast(locality~year)

samps_EFFCOV %>%
  filter(locality %in%
           c("DE_Bro","DE_Mun","FI_Aka","PA_li","TR_Yes","UA_Ode", "UA_od","VA_ch","WI_cp")) ->
  filtered_samps_for_analysis

### Run up to here to begin
#selected_pops = list(
#  all= unique(annot_samps$locality),
#  Bro = "DE_Bro",
#  Mun = "DE_Mun",
#  #Gim = "ES_Gim",
#  Aka = "FI_Aka",
#  li = "PA_li",
#  Yes = "TR_Yes",
#  #Kyi= "UA_Kyi",
#  Ode = c("UA_Ode", "UA_od"),
#  ch = "VA_ch",
#  cp = "WI_cp"
#)

#PCA_obj_ind_analysis = foreach(i=1:length(selected_pops), .combine = "rbind" )%do%{

  #message(names(selected_pops)[i])
  
  filtered_samps_for_analysis %>%
  dplyr::select(sampleId) %>% .$sampleId ->
  selected_samps

dat_for_Analysis %>%
  as.data.frame() %>%
  filter(rownames(.) %in% selected_samps) %>% 
  PCA(scale.unit = F, graph = F, ncp = 5) ->
  PCA_object

save(PCA_object, 
     file = "/scratch/yey2sn/Overwintering_ms/2.Temporal_Spatial_structure/RAW_PCA_obj_ind_analysis.allnalyses.Rdata")


  PCA_object$ind$coord %>%
    as.data.frame() %>%
    mutate(sampleId = rownames(.)
           #,
           #analysis_set = names(selected_pops)[i] 
           ) %>%
    left_join(filtered_samps_for_analysis) ->
    PCA_obj_ind_analysis
  
#} ## close loop

save(PCA_obj_ind_analysis, 
     file = "/scratch/yey2sn/Overwintering_ms/2.Temporal_Spatial_structure/PCA_obj_ind_analysis.allnalyses.Rdata")

load("/scratch/yey2sn/Overwintering_ms/2.Temporal_Spatial_structure/PCA_obj_ind_analysis.allnalyses.Rdata")

#save(PCA_object, 
#     file = "/scratch/yey2sn/Overwintering_ms/2.Temporal_Spatial_structure/PCA.object.all.Rdata")

###
###
###
###
###
###Build PCA on a chromosome by chromosome basis for other DEST samples.
### Run up to here to begin
selected_pops.chr = list(
  Cha = "VA_ch",
  Bro = "DE_Bro",
  Mun = "DE_Mun",
  #Gim = "ES_Gim",
  Aka = "FI_Aka",
  li = "PA_li",
  Yes = "TR_Yes",
  #Kyi= "UA_Kyi",
  Ode = c("UA_Ode", "UA_od"),
  Cp = "WI_cp"
)

### This object contains extra samps for analysis -- a more detailed metadata file
### 
chrs = c("2L","2R","3L", "3R")

chr.pcas.outer = foreach(k=1:length(chrs), .combine = "rbind" )%:% ## open outer nested loop 
  foreach(i=1:length(selected_pops.chr), .combine = "rbind")%do%{ ## open inner loop
    
    colnames(dat_for_Analysis) -> all_snps
    snps_chr = all_snps[grep(chrs[k], all_snps)]
    
    message(paste(names(selected_pops.chr)[i],chrs[k], sep = "-in-" ))
    
    filtered_samps_for_analysis %>%
      filter(locality %in% selected_pops.chr[[i]]) %>%
      dplyr::select(sampleId) %>% .$sampleId ->
      selected_samps
    
    message(selected_samps)
    
    ### subsample to pop of interest
    dat_for_Analysis %>%
      as.data.frame() %>% 
      filter(rownames(.) %in% selected_samps) %>% 
      .[,which(colnames(.) %in% snps_chr)] %>% 
      PCA(scale.unit = F, graph = F, ncp = 2) ->
      PCA_object
    
    PCA_object$ind$coord %>%
      as.data.frame() %>%
      mutate(sampleId = rownames(.),
             chr=chrs[k],
             analysis_set = names(selected_pops.chr)[i] ) %>%
      left_join(filtered_samps_for_analysis) -> tmp_obj
    tmp_obj
  } ## close inner loop

##  outer loop is closed by default using the %:% operator
##  outer loop is closed by default using the %:% operator

chr.pcas.outer$season[is.na(chr.pcas.outer$season)] = "time_series"

chr.pcas.outer %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    fill=year,
    #label=year,
    shape =season
  )) +
  #geom_text_repel(size = 1.5, max.overlaps = 10)+
  geom_point(size = 2) +
  scale_shape_manual(values = c(21,22,23,24)) +
  scale_fill_gradientn(
    colors=c("springgreen","cyan","blue","gold","red"),
    values=rescale(c(2011,2013,2015,2016,2018))
  ) +
  theme_bw() +
  facet_grid(city~chr) ->
  multipop_plot.chr
ggsave(multipop_plot.chr, file = "multipop_plot.chr.pdf", h = 8, w = 6.5)

save(chr.pcas.outer,
     file = "PCA_allpops_dims12.Rdata")


