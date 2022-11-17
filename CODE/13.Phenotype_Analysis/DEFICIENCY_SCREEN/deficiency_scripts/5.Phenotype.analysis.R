#startle response decay modeling
library(data.table)
library(tidyverse)
library(foreach)
library(lme4)
library(doParallel)
library(emmeans)
#library(ggrepel)
registerDoParallel(5)

setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Aug_2022_objects/")
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Oct_2022_objects/")

#goal- take our 90 sec aggregated data, and scale to activity per minute
#also scale b.a. data to activity per minute
#then take our window activity and mamke it scaled by - obs.act / mean(basal.act)

data2 = readRDS("data.smoothed90window30step.RDS")
data = readRDS("90win.30step")


flyid.vector = unique(data$fly.id)
#scale window.act and basal act to act/min
data$obs.act.min = data$window.activity *(2/3)
data$basal.act.min = data$basal_act / 60

#remove unneeded columns

data.zed2 = data %>% 
  group_by(fly.id, minute) %>% 
  mutate(zed.obs = (obs.act.min) / basal.act.min ) %>% 
  as.data.table(.)

#let's continue by finding the highest point of activity for each fly post stimulus, and having that be the initial time point for each fly.

highacttime = data.zed2 %>% 
  group_by(fly.id) %>% 
  filter(minute >= 0) %>% 
  mutate(highest.act = max(zed.obs, na.rm = T)) %>% #this is the highest activity
  filter(zed.obs == unique(highest.act)) %>% #this filters data down to each time point with highest activity
  filter(minute == min(minute)) %>% #now find first time point with the max
  summarise(high.act.timepoint = unique(minute)) #save this as a new variable
#merge this on
data.zed3 = merge(data.zed2, highacttime, by = "fly.id")
#now cut off previous time points
data.zed.cut = data.zed3 %>% 
  group_by(fly.id) %>% 
  filter(minute >= high.act.timepoint) %>% 
  as.data.table(.)
#next we want to change our time variable- instead of minutes since time zero, we want minutes since highest act(the start of our curve)

data.zed.cut = data.zed.cut %>% 
  group_by(fly.id, minute) %>% 
  mutate(scaled.minute = minute - high.act.timepoint) %>% 
  as.data.table(.)


flylongstartle = data.zed.cut %>% 
  group_by(fly.id) %>% 
  filter(zed.obs <= 1) %>% #filter otu to times in which activity returns to basal
  filter(scaled.minute == min(scaled.minute)) %>% # go to the first time point
  as.data.table(.)
#about 200 flys are sampled out- as ones that don't reach basal. 
#save this as the id and the phenotype (scaled minute)
savedata = flylongstartle %>% 
  mutate(sr.length = scaled.minute) %>% 
  dplyr::select(fly.id, sr.length) %>% 
  as.data.table(.)

#merge.startle duration on
phenotype.dt = merge(data.zed.cut, savedata, by = "fly.id" )
saveRDS(phenotype.dt, "startle.duration.11.9")
########################################################
###begin modeling- using mixed models###
########################################################
#split out metadata on genotype

phenotype.dt[,c("knockout", "other"):= tstrsplit(phenotype.dt$fly.id, "-")]

phenotype.dt[,c("DGRP", "f1.geno"):= tstrsplit(phenotype.dt$other, "\\.")]
phenotype.dt[,c("f1.background", "day", "week"):= tstrsplit(phenotype.dt$f1.geno, "_")]
phenotype.dt$inversion.st = ifelse(grepl("i", phenotype.dt$DGRP),"inverted", "standard")
phenotype.dt$f1.background = ifelse(grepl("B", phenotype.dt$f1.background),"balancer", "deficiency")
#create variable for msp status
status = data.frame (
  DGRP = c('di1', 'di2','di3','di4','di5','di6','ds1','ds2','ds3','ds4'),
  msp.status = c('alt','alt','ref','alt','alt','alt', 'ref','ref','ref','ref'),
  win5.bgd = c("derived", "derived", "ancestral", "ancestral", 'derived', "derived","ancestral","ancestral","ancestral","ancestral" ),
  win9.bgd = c("ancestral", 'ancestral', "derived", "derived", "ancestral", "derived", "ancestral","ancestral","ancestral","ancestral")
)

meltdata = merge(phenotype.dt, status, by = "DGRP")
saveRDS(phenotype.dt, "startle.data.RDS")
#meltdata = readRDS("startle.modeling.RDS")
#create geno label that uses ancestral or derived (type 1 or 2) and inversion status
data.2 = meltdata %>% 
  filter(knockout %in% c("k2b","k2a"))
data.5 = meltdata %>% 
  filter(knockout %in% c("k5a", "k5b"))
data.9 = meltdata %>% 
  filter(knockout == "k9a")
#create geno column that adds together inversion and ancestral derived
data.2 = data.2 %>% 
  mutate(geno = case_when(inversion.st == "inverted" ~ "inverted.derived",
                          inversion.st == "standard" ~ "standard.ancestral"))
data.5$geno = paste(data.5$inversion.st, data.5$win5.bgd, sep = ".")
data.9$geno = paste(data.9$inversion.st, data.9$win9.bgd, sep = ".")
meltdata = rbind(data.2, data.5, data.9) 

#remove time points before highpoint
data.merge = meltdata %>% 
  group_by(fly.id) %>% 
  filter(scaled.minute <= sr.length) %>% 
  mutate(geno = paste(geno, f1.background, sep = "-")) %>% 
  as.data.table(.)

#use loop to run emmeans to find slopes for each knockout
knockouts = c("k2a", "k2b", "k5a", "k5b", "k9a")
out = foreach(f = knockouts, .combine = "rbind") %do% {
 # f = knockouts[1]
  
  #make pairwise ~ (categorical factor), var = continues covariate
  model.em = glmer(formula = (zed.obs ~ scaled.minute * geno  + (1| fly.id) + (1|week)), data.merge[knockout == f])
  
  x = emtrends(model.em, pairwise ~ geno, var = "scaled.minute")#pbkrtest.limit = 6697
  #look for emmeans measure of intercepts. 
  
  #emmip(model.em, geno ~ scaled.minute, cov.reduce = range)
  x = x[[1]]
  em.data = as.data.table(x)
  em.data$f1.background = ifelse(grepl("balancer", em.data$geno),"balancer", "deficiency")
  em.data[,newgeno := tstrsplit(geno, "-")[[1]]]
  em.data$knockout = f
  em.data
}
saveRDS(out, "emmeans.data.RDS")

