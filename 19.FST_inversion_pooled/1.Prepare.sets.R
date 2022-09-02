##### Step 1. Prepare sets
##### 

rm(list = ls())

library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(magrittr)
library(doMC)
registerDoMC(4)
library(SeqArray)
library(lubridate)
library(tidyverse)
library(gtools)
library(poolfstat)
library(viridis)

##################
##################

### Create function for pair-wise comparions

create.comparisons.vector = function(samps.df, cont.c.lab, type.test){
  
  #samps.df = samps.ec.filtd.lab
  #cont.c.lab = "Cville.TYS"
  #type.test = "time"
  
  #cont.c.lab ="NoA.E"
  #type.test = "space"
  
  
  if(type.test == "time"){
    samps.df %>%
      filter(Continental_clusters == cont.c.lab) ->
      samps.tmp
  }

  if(type.test == "space"){
    samps.df %>%
      filter(Continental_clusters == cont.c.lab) %>%
      group_by(locality) %>%
      slice_head() ->
      samps.tmp
  }
  
  L = dim(samps.tmp)[1]
  
  comp_vector = combinations(
    L,
    2, 
    v=1:L,
    set=TRUE, 
    repeats.allowed=FALSE)
  
  print("Create combination vector")
  
  comp_vector %<>%
    as.data.frame() %>%
    mutate(day_diff = NA,
           year_diff = NA,
           same.locale = NA,
           same.season = NA
           )
  
  ##calculate day differences
  print("Loop to esrtimate time difference")
  
  for(i in 1:dim(comp_vector)[1]) {
    
    date1=samps.tmp$Date[comp_vector[i,1]]
    date2=samps.tmp$Date[comp_vector[i,2]]
    comp_vector$day_diff[i] = abs(as.numeric(date1-date2))
    
    year1=samps.tmp$year[comp_vector[i,1]]
    year2=samps.tmp$year[comp_vector[i,2]]
    comp_vector$year_diff[i] = abs(as.numeric(year1-year2))
    
    season1=samps.tmp$season[comp_vector[i,1]]
    season2=samps.tmp$season[comp_vector[i,2]]
    if(is.na(season1) | is.na(season2) ){
    comp_vector$same.season[i] = NA} else{
      comp_vector$same.season[i] = ifelse(season1 == season2, 
                                          "yes", "no")}
  
    comp_vector$pop1[i] = samps.tmp$city[comp_vector[i,1]]
    comp_vector$pop2[i] = samps.tmp$city[comp_vector[i,2]]
    
    comp_vector$samp1[i] = samps.tmp$sampleId[comp_vector[i,1]]
    comp_vector$samp2[i] = samps.tmp$sampleId[comp_vector[i,2]]
    
    comp_vector$same.locale[i] = ifelse(comp_vector$pop1[i] == comp_vector$pop2[i], 
                                        "yes", "no")
  
  }
  
  return(comp_vector)
  
}

##################
##################

samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
samps[nchar(month)==1, month:=paste("0", month, sep="")]
samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
samps[nchar(day)==1, day:=paste("0", day, sep="")]
samps <- samps[set!="dgn"]
samps[,Date:=date(paste(year, month, day, sep="-"))]

inmeta="/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST_EC_metadata.Rdata"
load(inmeta)
samps_EFFCOV %>% 
  dplyr::select(sampleId, MeanEC) ->
  sampsEC

file.clusters = "/scratch/yey2sn/Overwintering_ms/19.inv.fst/DEST_Sample_clusters.txt"
clust.asg = fread(file.clusters)

clust.asg %>%
  dplyr::select(sampleId, Continental_clusters) ->
  samps.clust

###
samps %>%
  left_join(sampsEC) %>%
  left_join(samps.clust) %>%
  filter(MeanEC > 28) ->
  samps.ec.filtd
###

samps.ec.filtd %>%
  mutate(Continental_clusters = case_when( 
   locality %in% c("VA_ch", "FL_ho", "GA_at", "GA_ha", "MA_la", "ME_bo", "NC_ra", "NY_it", "PA_li", "PA_st", "SC_eu") ~ "NoA.E",
   TRUE ~ Continental_clusters
    ),
   Continental_clusters = case_when( 
     set %in% c("CvilleSet") ~ "Cville.TYS",
     TRUE ~ Continental_clusters
   )) ->
  samps.ec.filtd.lab

### --- preapre core 20 samps
core20.samps <- read.csv("./core20_samps.csv", row.names = F)
core20.samps %<>%
  mutate(sampleId = rownames(.))

samps.ec.filtd.lab %>%
  filter(sampleId %in% core20.samps$sampleId ) %>%
  mutate(Continental_clusters = "Core20") ->
  core20.df

### 
samps.ec.filtd.lab %>% .$Continental_clusters %>% table

####
create.comparisons.vector(samps.ec.filtd.lab, "Cville.TYS", "time") %>% 
  filter(year_diff %in% 0:1) %>%
  mutate(comp.set = "time.cville") ->
  time.cville.comps

create.comparisons.vector(samps.ec.filtd.lab, "NoA.E", "space") %>% 
  filter(same.locale == "no")  %>%
  mutate(comp.set = "space.NoA.E") ->
  space.NoA.E.comps

create.comparisons.vector(samps.ec.filtd.lab, "1.Europe_W", "time") %>% 
  filter(year_diff %in% 0:1) %>% 
  filter(same.locale == "yes")  %>%
  mutate(comp.set = "time.eu.w") ->
  time.eu.w.comps

create.comparisons.vector(samps.ec.filtd.lab, "3.Europe_E", "time") %>% 
  filter(year_diff %in% 0:1) %>% 
  filter(same.locale == "yes")  %>%
  mutate(comp.set = "time.eu.e") ->
  time.eu.e.comps

create.comparisons.vector(samps.ec.filtd.lab, "1.Europe_W", "space") %>% 
  filter(same.locale == "no")  %>%
  mutate(comp.set = "space.eu.w") ->
  space.eu.w.comps

create.comparisons.vector(samps.ec.filtd.lab, "3.Europe_E", "space") %>% 
  filter(same.locale == "no")  %>%
  mutate(comp.set = "space.eu.e") ->
  space.eu.e.comps


create.comparisons.vector(core20.df, "Core20", "time") %>% 
  filter(same.locale == "yes")  %>%
  filter(same.season == "no")  %>% 
  mutate(comp.set = "core20.seasonality") ->
  core20.seasonality.comps

master.comp.list =
rbind(time.cville.comps,
space.NoA.E.comps,
time.eu.w.comps,
time.eu.e.comps,
space.eu.w.comps,
space.eu.e.comps,
core20.seasonality.comps)

write.table(master.comp.list, 
            file = "master.comp.list.fst.txt", 
            append = FALSE, 
            quote = F, sep = "\t",
            eol = "\n", na = "NA", 
            dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

