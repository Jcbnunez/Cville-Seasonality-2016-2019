#SNP-wise analysis
# Verison 2 - August 19, 2022
# 
#
#
## plot analysis
## 


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

load("cville.V2.snp.wise.df.Rdata")
##snp.wise.tmp.f.df

#####
load("/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt_25percMissing.Rdata")
head(snp.dt)
#####
#####

#### import genomic data
genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)

samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
samps[nchar(month)==1, month:=paste("0", month, sep="")]
samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
samps[nchar(day)==1, day:=paste("0", day, sep="")]
samps <- samps[set!="dgn"]
samps[,Date:=date(paste(year, month, day, sep="-"))]

samps %>%
  dplyr::select(samp1 = sampleId, locality_1 = locality,  year_1 = year, Month_1 = month, season_1 = season ) ->
  samp1

samps %>%
  dplyr::select(samp2 = sampleId, locality_2 = locality,  year_2 = year,  Month_2 = month,  season_2 = season ) ->
  samp2

### get subsample of data to work on
seqResetFilter(genofile)
#seqSetFilter(genofile, sample.id=samps.cville$sampleId)
seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snp.dt$id)
snp.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

###

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
  filter(rnp < 0.05) %>%
  mutate(win = case_when(
    pos > 2800000 & pos < 3200000 ~ "win_3.1",
    pos > 4470000 & pos < 4870000 ~ "win_4.7",
    pos > 4920000 & pos < 5320000 ~ "win_5.1",
    pos > 6000000 & pos < 6400000 ~ "win_6.1",
    pos > 6600000 & pos < 7000000 ~ "win_6.8",
    pos > 9400000 & pos < 9800000 ~ "win_9.6"
  )) %>%
  mutate(glm_snp = paste(chr, pos , sep = "_")) %>%
  filter(!is.na(win)) ->
  glm.out.cville.win

glm.out.cville.win %<>%
  left_join(snp.dt[,c("chr", "pos", "cm_mb", "af")]) %>%
  mutate(snp = paste(chr, pos, sep = "_"))

##glm.outliers
glm.out.cville.win %>%
  left_join(snp.wise.tmp.f.df) %>%
  group_by(win) %>%
  slice_min(p_lrt) -> glm.out.cville.win.annot

glm.out.cville.win.annot %>%
  group_by(snp) %>%
  slice_head(n=1) %>%
  .[complete.cases(.$pan.af),] %>%
  dplyr::select(snp, pan.af, cm_mb) -> glms.ids


##########
##########
##########
##########
##make pool of controls
glm.out.cville %>%
  filter(chr != "2L") %>%
  filter(invName != "none") %>%
  filter(perm == 0) %>%
  filter(rnp > 0.5) %>%
  left_join(snp.dt[,c("chr", "pos", "cm_mb", "af")]) %>%
  mutate(snp = paste(chr, pos, sep = "_")) ->
  pool_of_controls

pool_of_controls %>%
  left_join(snp.wise.tmp.f.df) %>% filter(!is.na(snp.FST)) ->
  pool_of_controls.annot

pool_of_controls.annot %>%
  group_by(snp) %>%
  slice_head(n=1) %>%
  .[complete.cases(.$pan.af),] %>%
  dplyr::select(snp, pan.af, cm_mb) -> controls.ids

#####
#####
#####
in.df=glms.ids
cont.id=controls.ids
chosen_controls.cville = foreach(i = 1:dim(in.df)[1],
                                 .combine = "rbind",
                                 .errorhandling ="remove")%do%{
                                   
                                   message(i)
                                   in.df[i,] ->
                                     tmp.achor
                                   
                                   cont.id %>%
                                     filter(cm_mb > tmp.achor$cm_mb-0.05 & cm_mb < tmp.achor$cm_mb+0.05,
                                            pan.af > tmp.achor$pan.af-0.01 & pan.af < tmp.achor$pan.af+0.01
                                     ) %>%
                                     group_by(snp) %>%
                                     .[10,]  ->
                                     select.control
                                   
                                   data.frame(tmp.achor, control_snp =  select.control$snp )
                                   
                                 }

chosen_controls.cville %<>%
  filter(!is.na(control_snp)) 
names(chosen_controls.cville)[1] = "glm_snp"

#### modify fst files
glm.out.cville.win.annot %>%
  dplyr::select(chr, pos, win, glm_snp, glm.fst = snp.FST, samp1, samp2, set) -> glm.to.merge

pool_of_controls.annot %>%
  dplyr::select(control_snp = snp, control.fst = snp.FST, samp1, samp2, set) -> control.to.merge

###
chosen_controls.cville %>%
  left_join(glm.to.merge) %>% 
  left_join(control.to.merge) %>% 
  left_join(samp1) %>%
  left_join(samp2) ->
  glms.and.control.matched.fsts

save(glms.and.control.matched.fsts, file = "glms.and.control.matched.fsts.Rdata")
load("glms.and.control.matched.fsts.Rdata")

glms.and.control.matched.fsts %>%
  filter(year_1 %in% 2016:2018 & year_2 %in% 2016:2018 ) %>%
  mutate(dy = abs(year_1-year_2),
         dm = abs(as.numeric(Month_1)-as.numeric(Month_2))) %>% 
  filter(dy %in% 0:1) ->
  dat.for.plot

####
dat.for.plot %>%
  group_by( dy, glm_snp, win) %>%
  summarize(mean.fst.glm = mean(glm.fst, na.rm = T),
            sd.fst.glm = sd(glm.fst, na.rm = T),
            mean.fst.control = mean(control.fst, na.rm = T),
            sd.fst.control = sd(control.fst, na.rm = T),
  ) -> dat.for.plot.sum

dat.for.plot.sum %>%
  ggplot() +
  geom_violin(aes(x=as.factor(dy),
                    y=glm.fst)) +
  geom_point(aes(x = as.factor(dy),
                  y= control.fst, 
                 color = dm)) +
  facet_grid(dm~win)->
  density.plot

ggsave(density.plot, file ="density.plot.png", h = 4, w = 7)


####
dat.for.plot %>%
  group_by( dy, dm, win) %>%
  summarize(mean.fst.glm = median(glm.fst, na.rm = T),
            mean.fst.control = median(control.fst, na.rm = T),
            ) %>% 
  melt(id = c(#"glm_snp", 
              "dy", "dm", "win")) %>% 
  ggplot(aes(
    #shape = variable,
    #linetype = variable,
    color= as.factor(dy),
    x= as.factor(dy) ,
    y= value
  )) +
  geom_violin() +
  ylab("mean FST") +
  xlab("month difference") +
  geom_hline(yintercept = 0) +
  geom_point() +
    facet_grid(dm~win)->
    fst.matched.vio
  
  ggsave(fst.matched.vio, file ="fst.matched.vio.png", h = 4, w = 7)
  