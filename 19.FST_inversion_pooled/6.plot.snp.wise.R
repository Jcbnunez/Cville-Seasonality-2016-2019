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

## samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
## samps[,c("sampleId", "locality", "season", "continent", "long", "lat", "year", "yday")] -> samps.tmp
## samp1 = samps.tmp
## names(samp1) = c("samp1", "locality_1", "season_1", "continent_1", "long_1", "lat_1", "year_1", "yday_1")
## samp2 = samps.tmp
## names(samp2) = c("samp2", "locality_2", "season_2", "continent_2", "long_2", "lat_2", "year_2", "yday_2")
####
####
extra_data <- "/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_2L.ECfiltered.Rdata"
load(extra_data)
head(filtered_samps_for_analysis)

filtered_samps_for_analysis %>%
  dplyr::select(samp1 = sampleId, locality_1 = locality,  year_1 = year, NC_1 = MeanEC, Month_1 = Month, season_1 = season ) ->
  samp1

filtered_samps_for_analysis %>%
  dplyr::select(samp2 = sampleId, locality_2 = locality,  year_2 = year, NC_2 = MeanEC, Month_2 = Month,  season_2 = season ) ->
  samp2


##### CREATE FUNCTION
##### ##### CREATE FUNCTION
##### CREATE FUNCTION

estimate.fs.means = function(in.df, LOC.tar){
  in.df  %>%
    mutate(het.bin.round = round(w.Het.snp, 1)) %>% 
    mutate(dy = abs(year_1-year_2)) %>%
    mutate(dm = abs(Month_1-Month_2)) %>%
    filter(dy %in% 0:1) %>%
    filter(locality_1 == locality_2 ) %>%
    filter(locality_1 %in% LOC.tar) ->
    markers.het.fst.cville.select
  
  if(LOC.tar == "VA_ch"){
    markers.het.fst.cville.select %<>%
      filter(year_1 %in% 2016:2018 & year_2 %in% 2016:2018 )
  }
  
  markers.het.fst.cville.select %>% 
    mutate(rnp.sig = case_when(
      p_lrt < 0.0000001 ~ "1e-07",
      p_lrt < 0.000001 ~ "1e-06",
      p_lrt < 0.00001 ~ "1e-05",
      p_lrt < 0.0001 ~ "1e-04",
      p_lrt < 0.001 ~ "1e-03",
      p_lrt < 0.01 ~ "1e-02",
      p_lrt < 0.05 ~ "5e-02",
      #TRUE ~ "Controls"
    ) )->
    tmp
  
  tmp$rnp.sig[is.na(tmp$rnp.sig)] = "Controls"
  
  tmp$rnp.sig = factor(tmp$rnp.sig, levels = c("Controls", "5e-02", "1e-02", "1e-03", "1e-04", "1e-05",
                                               "1e-06","1e-07") )
  
  return(tmp)
}

##### CREATE FUNCTION



####
####
base <- "/project/berglandlab/alan/environmental_ombibus_global"
#########
model = "temp.max;2;5.Cville"
file.cvile <- paste(base, model, 
                    paste(model, ".glmRNP.Rdata", sep = ""), 
                    sep = "/" )
print(file.cvile)
glm.out.cville <- get(load(file.cvile))
glm.out.cville %>%
  filter(chr == "2L") %>%
  filter(perm == 0) %>%
  mutate(snp = paste(chr, pos , sep = "_")) %>%
  dplyr::select(snp, chr, pos, rnp, p_lrt) ->
  glm.out.cville.rnp

## load matched controls
load("./glm.and.matchedControls.Rdata")
matched.inv.cville <- macthed.controls.Inv %>% mutate(snp = paste(chr, pos, sep = "_") )
macthed.no.inv.cville <- macthed.controls.noInv %>% mutate(snp = paste(chr, pos, sep = "_") )

##file.name <- "/scratch/yey2sn/Overwintering_ms/19.inv.fst/glm.snps.cville.snp.wise.tmp.f.df.Rdata"

file.names <- system(" ls | grep -e 'snp.wise' | grep 'cville' ", intern = T)
#=1
markers.het.fst.cville = foreach(i = 1:length(file.names), 
                          .combine = "rbind")%do%{
####
                            file.name = file.names[i]
                            snp.wise.name = tstrsplit(file.name, "\\.")
                            pop = ifelse(snp.wise.name[[1]] == "glm", snp.wise.name[[3]], snp.wise.name[[4]]  )
                            type.of.marker = ifelse(snp.wise.name[[1]] == "glm", "glm", snp.wise.name[[3]]  )
                            
                            message(paste(i, file.name, pop, type.of.marker, sep = " | " ))
                            
load(file.name)
snp.wise.tmp.f.df %>%
  left_join(samp1) %>%
  left_join(samp2) %>% 
  filter(NC_1 > 10 & NC_2 > 10) %>% 
  filter(!is.na(snp.FST)) %>%
  left_join(glm.out.cville.rnp) -> snp.wise.tmp.f.df.annot

snp.wise.tmp.f.df.annot %>%
  mutate(pop = pop, type.of.marker = type.of.marker)  ->
  snp.wise.tmp.f.ag.annot

snp.wise.tmp.f.ag.annot %<>%
  mutate(win = case_when(
    pos > 2800000 & pos < 3200000 ~ "win_3.1",
    pos > 4470000 & pos < 4870000 ~ "win_4.7",
    pos > 4920000 & pos < 5320000 ~ "win_5.1",
    pos > 6000000 & pos < 6400000 ~ "win_6.1",
    pos > 6600000 & pos < 7000000 ~ "win_6.8",
    pos > 9400000 & pos < 9800000 ~ "win_9.6"
  ))

return(snp.wise.tmp.f.ag.annot)
}

markers.het.fst.cville %>% 
  filter(!is.na(win)) -> fst.glm.wins

fst.glm.wins %>%
  left_join(matched.inv.cville) %>% head


markers.het.fst.cville %>% 
  filter(set == "macthed.controls.noInv") %>%
  left_join(macthed.no.inv.cville) %>% head
  filter(!is.na(matched.to)) %>% head
  







####
markers.het.fst.cville.select =
estimate.fs.means( markers.het.fst.cville ,"VA_ch")

### plot VA
markers.het.fst.cville.select %>%
  filter(Month_1 %in% c(8,11) &  Month_2 %in% c(8,11)) %>% 
  group_by(dy, dm, rnp.sig) %>%
  summarize(me.fst = mean(snp.FST),
            sd.fst = sd(snp.FST)) %>%
  ggplot(
    aes( y=(me.fst),
         ymin=me.fst-sd.fst,
         ymax=me.fst+sd.fst,
         x = as.factor(rnp.sig),
         color = as.factor(dy), 
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(alpha = 0.5,
                position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) +
  coord_flip() +
  xlab("GLM P treshold") +
  ylab("Mean FST +/- SD") +
  ggtitle("August vs November | within vs Overwintering", subtitle = "Cville") + 
  theme_classic() +
  facet_grid(paste(dm) ~ .) ->
  fst.density2.cville

ggsave(fst.density2.cville, file = "fst.density2.cville.pdf", h = 4, w = 5)

########
markers.het.fst.PAli =
  estimate.fs.means( markers.het.fst.cville ,"PA_li")

### plot PA
markers.het.fst.PAli %>%
  filter(season_1 != season_2) %>% 
  group_by(dy, rnp.sig) %>%
  summarize(me.fst = mean(snp.FST),
            sd.fst = sd(snp.FST)) %>%
  ggplot(
    aes( y=(me.fst),
         ymin=me.fst-sd.fst,
         ymax=me.fst+sd.fst,
         x = as.factor(rnp.sig),
         color = as.factor(dy), 
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(alpha = 0.5,
                position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) +
  coord_flip() +
  xlab("GLM P treshold") +
  ylab("Mean FST +/- SD") +
  ggtitle("Spring vs Fall | within vs Overwintering", subtitle = "Lindvilla") + 
  theme_classic()  ->
  fst.density2.PA

ggsave(fst.density2.PA, file = "fst.density2.PA.pdf", h = 4, w = 5)

#########
#########
#########
#########

######
base <- "/project/berglandlab/alan/environmental_ombibus_global"

model = "temp.ave;9;3.Europe_E"
file.eue <- paste(base, model, 
                    paste(model, ".glmRNP.Rdata", sep = ""), 
                    sep = "/" )
print(file.eue)
glm.out.eue <- get(load(file.eue))
glm.out.eue %>%
  filter(chr == "2L") %>%
  filter(perm == 0) %>%
  mutate(snp = paste(chr, pos , sep = "_")) %>%
  dplyr::select(snp, rnp, p_lrt) ->
  glm.out.eue.rnp

file.names <- system(" ls | grep -e 'snp.wise' | grep 'EUE' ", intern = T)
#=1
markers.het.fst.eue = foreach(i = 1:length(file.names), 
                                 .combine = "rbind")%do%{
                                   ####
                                   file.name = file.names[i]
                                   snp.wise.name = tstrsplit(file.name, "\\.")
                                   pop = ifelse(snp.wise.name[[1]] == "glm", snp.wise.name[[3]], snp.wise.name[[4]]  )
                                   type.of.marker = ifelse(snp.wise.name[[1]] == "glm", "glm", snp.wise.name[[3]]  )
                                   
                                   message(paste(i, file.name, pop, type.of.marker, sep = " | " ))
                                   
                                   load(file.name)
                                   snp.wise.tmp.f.df %>%
                                     left_join(samp1) %>%
                                     left_join(samp2) %>% 
                                     filter(NC_1 > 28 & NC_2 > 28) %>% 
                                     filter(!is.na(snp.FST)) %>%
                                     left_join(glm.out.eue.rnp) -> snp.wise.tmp.f.df.annot
                                   
                                   snp.wise.tmp.f.df.annot %>%
                                     mutate(pop = pop, type.of.marker = type.of.marker)  ->
                                     snp.wise.tmp.f.ag.annot
                                   
                                   return(snp.wise.tmp.f.ag.annot)
                                 }

###
markers.het.fst.Odesa =
  estimate.fs.means( markers.het.fst.eue ,"UA_Ode")

### plot Ode
markers.het.fst.Odesa %>%
  filter(Month_1 %in% c(7, 10) & Month_2 %in% c(7, 10) ) %>% 
  group_by(dy, rnp.sig, dm) %>%
  summarize(me.fst = mean(snp.FST),
            sd.fst = sd(snp.FST)) %>% 
  ggplot(
    aes( y=(me.fst),
         ymin=me.fst-sd.fst,
         ymax=me.fst+sd.fst,
         x = as.factor(rnp.sig),
         color = as.factor(dy), 
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(alpha = 0.5,
                position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) +
  coord_flip() +
  xlab("GLM P treshold") +
  ylab("Mean FST +/- SD") +
  ggtitle("Spring vs Fall | within vs Overwintering", subtitle = "Odessa") + 
  facet_grid(dm ~ .) +
  theme_classic()  ->
  fst.density2.Ode

ggsave(fst.density2.Ode, file = "fst.density2.Ode.pdf", h = 4, w = 5)

#########
#########
#########
#########

######
model = "humidity.ave;8;1.Europe_W"
file.euw <- paste(base, model, 
                  paste(model, ".glmRNP.Rdata", sep = ""), 
                  sep = "/" )
print(file.euw)
glm.out.euw <- get(load(file.euw))
glm.out.euw %>%
  filter(chr == "2L") %>%
  filter(perm == 0) %>%
  mutate(snp = paste(chr, pos , sep = "_")) %>%
  dplyr::select(snp, rnp, p_lrt) ->
  glm.out.euw.rnp


####
file.names <- system(" ls | grep -e 'snp.wise' | grep 'EUW' ", intern = T)
markers.het.fst.euw = foreach(i = 1:length(file.names), 
                              .combine = "rbind")%do%{
                                ####
                                file.name = file.names[i]
                                snp.wise.name = tstrsplit(file.name, "\\.")
                                pop = ifelse(snp.wise.name[[1]] == "glm", snp.wise.name[[3]], snp.wise.name[[4]]  )
                                type.of.marker = ifelse(snp.wise.name[[1]] == "glm", "glm", snp.wise.name[[3]]  )
                                
                                message(paste(i, file.name, pop, type.of.marker, sep = " | " ))
                                
                                load(file.name)
                                snp.wise.tmp.f.df %>%
                                  left_join(samp1) %>%
                                  left_join(samp2) %>% 
                                  filter(NC_1 > 28 & NC_2 > 28) %>% 
                                  filter(!is.na(snp.FST)) %>%
                                  left_join(glm.out.euw.rnp) -> snp.wise.tmp.f.df.annot
                                
                                snp.wise.tmp.f.df.annot %>%
                                  mutate(pop = pop, type.of.marker = type.of.marker)  ->
                                  snp.wise.tmp.f.ag.annot
                                
                                return(snp.wise.tmp.f.ag.annot)
                              }

###
markers.het.fst.euw$locality_1 %>% table
markers.het.fst.samp.euw =
  estimate.fs.means( markers.het.fst.euw ,"DE_Mun")

### plot Ode
markers.het.fst.samp.euw %>%
  #filter(Month_1 %in% c(6, 9) & Month_2 %in% c(6, 9) ) %>% 
  group_by(dy, rnp.sig, dm) %>%
  summarize(me.fst = mean(snp.FST),
            sd.fst = sd(snp.FST)) %>% 
  ggplot(
    aes( y=(me.fst),
         ymin=me.fst-sd.fst,
         ymax=me.fst+sd.fst,
         x = as.factor(rnp.sig),
         color = as.factor(dy), 
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(alpha = 0.5,
                position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) +
  coord_flip() +
  xlab("GLM P treshold") +
  ylab("Mean FST +/- SD") +
  ggtitle("Spring vs Fall | within vs Overwintering", subtitle = "test") + 
  facet_grid(dm ~ .) +
  theme_classic()  ->
  fst.density2.samp.euw

ggsave(fst.density2.samp.euw, file = "fst.density2.samp.euw.pdf", h = 4, w = 5)








LOC.tar = markers.het.fst.EUE$locality_1 %>% table %>% names

in.obj=markers.het.fst.EUE
eue.thsd.dat = foreach(LOC.tar = LOC.tar, .combine = "rbind", .errorhandling = "remove")%do%{
  
  message(LOC.tar)
  
  in.obj %>%
    mutate(het.bin.round = round(w.Het.snp, 1)) %>% 
    mutate(dy = abs(year_1-year_2)) %>%
    #filter(dy == 0) %>%
    filter(rd.sum > 120) %>%
    filter(locality_1 == locality_2 ) %>%
    filter(locality_1 == LOC.tar) ->
    tmp.dat.pre
  
  het.thsd.dat = foreach(het = seq(from = 0.1, to = 0.5, by = 0.1), 
                         .combine = "rbind",
                         .errorhandling = "remove")%do%{
                           
                           message(het)
                           tmp.dat.pre %>%
                             filter(het.bin.round == het) ->
                             tmp.dat
                           
                           message("tukey hsd")
                           
                           lm(snp.FST ~  set, data = tmp.dat) %>% aov %>% TukeyHSD(conf.level = 0.99) -> 
                             TukeyHSD.tmp.dat
                           
                           TukeyHSD.tmp.dat$set %>%
                             as.data.frame() %>% 
                             mutate(comp = rownames(.),
                                    pop = LOC.tar,
                                    het = het) -> tmp.out
                           
                           return(tmp.out)
                         }
  return(het.thsd.dat)
  
}

save(eue.thsd.dat, file = "eue.thsd.dat.Rdata")

###
######

file.names <- system(" ls | grep -e 'snp.wise' | grep 'EUW' ", intern = T)
#=1
markers.het.fst.EUW = foreach(i = 1:length(file.names), 
                              .combine = "rbind")%do%{
                                ####
                                file.name = file.names[i]
                                snp.wise.name = tstrsplit(file.name, "\\.")
                                pop = ifelse(snp.wise.name[[1]] == "glm", snp.wise.name[[3]], snp.wise.name[[4]]  )
                                type.of.marker = ifelse(snp.wise.name[[1]] == "glm", "glm", snp.wise.name[[3]]  )
                                
                                message(paste(i, file.name, pop, type.of.marker, sep = " | " ))
                                
                                load(file.name)
                                snp.wise.tmp.f.df %>%
                                  left_join(samp1) %>%
                                  left_join(samp2) -> snp.wise.tmp.f.df.annot
                                
                                snp.wise.tmp.f.df.annot %>%
                                  mutate(comp.type = case_when( locality_1 == "VA_ch" & locality_2 == "VA_ch" ~ "VA_comp",
                                                                TRUE ~ "notVaComp")) ->
                                  snp.wise.tmp.f.df.annot2
                                
                                ### generate sets
                                snp.wise.tmp.f.df.annot2 %>%
                                  .[complete.cases(.$snp.FST),] %>%
                                  group_by(comp.type, snp, set) %>%
                                  mutate(mean.Het = mean(w.Het.snp),
                                         mean.gAF = mean(pan.af),
                                         mean.snp.Q2 = mean(snp.Q2),
                                         mean.snp.fst = mean(snp.FST)
                                  ) -> 
                                  snp.wise.tmp.f.ag
                                
                                snp.wise.tmp.f.ag %>%
                                  mutate(pop = pop, type.of.marker = type.of.marker)  ->
                                  snp.wise.tmp.f.ag.annot
                                
                                return(snp.wise.tmp.f.ag.annot)
                              }

LOC.tar = markers.het.fst.EUW$locality_1 %>% table %>% names

in.obj=markers.het.fst.EUW
euw.thsd.dat = foreach(LOC.tar = LOC.tar, .combine = "rbind", .errorhandling = "remove")%do%{
  
  message(LOC.tar)
  
  in.obj %>%
    mutate(het.bin.round = round(w.Het.snp, 1)) %>% 
    mutate(dy = abs(year_1-year_2)) %>%
    #filter(dy == 0) %>%
    filter(rd.sum > 120) %>%
    filter(locality_1 == locality_2 ) %>%
    filter(locality_1 == LOC.tar) ->
    tmp.dat.pre
  
  het.thsd.dat = foreach(het = seq(from = 0.1, to = 0.5, by = 0.1), 
                         .combine = "rbind",
                         .errorhandling = "remove")%do%{
                           
                           message(het)
                           tmp.dat.pre %>%
                             filter(het.bin.round == het) ->
                             tmp.dat
                           
                           message("tukey hsd")
                           
                           lm(snp.FST ~  set, data = tmp.dat) %>% aov %>% TukeyHSD(conf.level = 0.99) -> 
                             TukeyHSD.tmp.dat
                           
                           TukeyHSD.tmp.dat$set %>%
                             as.data.frame() %>% 
                             mutate(comp = rownames(.),
                                    pop = LOC.tar,
                                    het = het) -> tmp.out
                           
                           return(tmp.out)
                         }
  return(het.thsd.dat)
  
}

save(euw.thsd.dat, file = "euw.thsd.dat.Rdata")
 ####
 ####
 ####

eue.thsd.dat <- get(load("eue.thsd.dat.Rdata"))
euw.thsd.dat <- get(load("euw.thsd.dat.Rdata"))

rbind(eue.thsd.dat, euw.thsd.dat) %>%
  filter(pop %in% c("UA_Ode", "DE_Bro", "DE_Mun", "FI_Aka", "TR_Yes")  ) %>%
  ggplot(aes(
    x=as.factor(het),
    y=diff,
    ymin=lwr,
    ymax=upr,
    color =comp
  )) +
  geom_errorbar(position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) +
  facet_grid(~pop) +
  theme(legend.pos = "bottom") +
  ylab("GLM higher fst <--> Controls higher fst") +
  geom_hline(yintercept = 0) ->
  plot.fst.thsd.eu

ggsave(plot.fst.thsd.eu, file = "plot.fst.thsd.eu.pdf", w = 8, h = 3 )
