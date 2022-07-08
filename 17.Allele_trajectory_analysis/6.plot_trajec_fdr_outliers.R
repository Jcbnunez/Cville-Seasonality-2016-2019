#### Plot trajectories of GLM outliers
#### 
#### 

library(tidyverse)
library(magrittr)
library(data.table)
library(car)
library(DescTools)
library(foreach)
library(doMC)
library(patchwork)
library(ggbeeswarm)
library(reshape2)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(viridis)
library(lubridate)
library(forcats)
library(viridis)
library(SeqArray)
library(tidyverse)
library(gmodels)
library(scatterpie)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(weathermetrics)
registerDoMC(2)

###
base <- "/project/berglandlab/alan/environmental_ombibus_global"

#########
model = "temp.max;2;5.Cville"

file.cvile <- paste(base, model, 
                    paste(model, ".glmRNP.Rdata", sep = ""), 
                    sep = "/" )

print(file.cvile)

glm.out <- get(load(file.cvile))

####

glm.out %>%
  filter(perm == 0) %>%
  mutate(fdr.p = p.adjust( p_lrt, "fdr")) ->
  glm.out.padjust

glm.out.padjust %<>%
  mutate(Win = case_when(
                pos > 2000000 & pos < 2400000 ~ "left",
                pos > 12900000 & pos < 13300000 ~ "right",
                pos > 2800000 & pos < 3200000 ~ "w3.1",
                pos > 4470000 & pos < 4870000 ~ "w4.7",
                pos > 4920000 & pos < 5320000 ~ "w5.1",
                pos > 6000000 & pos < 6400000 ~ "w6.1",
                pos > 6600000 & pos < 7000000 ~ "w6.8",
                pos > 9400000 & pos < 9800000 ~ "w9.6"))

glm.out.padjust %>%
  filter(fdr.p <= 0.1) %>%
  group_by(Win) %>%
  summarise(N = n())

glm.out.padjust %>%
  filter(!is.na(Win)) %>%
  filter(fdr.p <= 0.1) %>%
  .$pos -> tag.pos.cands

load("/scratch/yey2sn/Overwintering_ms/7.LD/merged.ld.Rdata")
### Add the bp distance between SNPs
ld_df %<>% 
  mutate(BP_diff = abs(BP_A-BP_B)) 

ld_df %>%
  filter(BP_diff != 0) %>%
  filter(BP_diff < 200000 & CHR_A == "2L" & CHR_B == "2L") %>%
  filter(BP_A %in% tag.pos.cands | BP_B %in% tag.pos.cands) %>%
  filter(R2 > 0.6) -> ld_df.tags.cands

ld_df.tags.cands %>%
  mutate(Win = case_when(
    BP_A > 2000000 & BP_A < 2400000 ~ "left",
    BP_A > 12900000 & BP_A < 13300000 ~ "right",
    BP_A > 2800000 & BP_A < 3200000 ~ "w3.1",
    BP_A > 4470000 & BP_A < 4870000 ~ "w4.7",
    BP_A > 4920000 & BP_A < 5320000 ~ "w5.1",
    BP_A > 6000000 & BP_A < 6400000 ~ "w6.1",
    BP_A > 6600000 & BP_A < 7000000 ~ "w6.8",
    BP_A > 9400000 & BP_A < 9800000 ~ "w9.6")) %>%
  filter(!is.na(Win)) %>%
  group_by(BP_A, Win) %>%
  summarise(N = n()) %>%
  group_by(Win) %>%
  slice_max(N, with_ties = F) -> ld_df.tags.cands.finalists

ld_df.tags.cands %>%
  filter(BP_A %in% ld_df.tags.cands.finalists$BP_A | BP_B %in% ld_df.tags.cands.finalists$BP_A) -> ld_df.tags.cands.finalists.all

unique(
ld_df.tags.cands.finalists.all$BP_A,
ld_df.tags.cands.finalists.all$BP_B) -> haplotags.snps

####


### load meta-data file
samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")

### open GDS for common SNPs (PoolSNP)
genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds")

### common SNP.dt
snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                     pos=seqGetData(genofile, "position"),
                     nAlleles=seqGetData(genofile, "$num_allele"),
                     id=seqGetData(genofile, "variant.id"))
snp.dt <- snp.dt[nAlleles==2]
seqSetFilter(genofile, snp.dt$id)

### function
getData <- function(chr="2L", start=14617051, end=14617051) {
  # chr="2L"; start=14617051; end=14617051
  
  ### filter to target
  snp.tmp <- data.table(chr=chr, pos=start:end)
  setkey(snp.tmp, chr, pos)
  setkey(snp.dt, chr, pos)
  seqSetFilter(genofile, variant.id=snp.dt[J(snp.tmp), nomatch=0]$id)
  
  ### get annotations
  #message("Annotations")
  tmp <- seqGetData(genofile, "annotation/info/ANN")
  len1 <- tmp$length
  len2 <- tmp$data
  
  snp.dt1 <- data.table(len=rep(len1, times=len1),
                        ann=len2,
                        id=rep(snp.dt[J(snp.tmp), nomatch=0]$id, times=len1))
  
  # Extract data between the 2nd and third | symbol
  snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
  snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]
  
  # Collapse additional annotations to original SNP vector length
  snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=","), gene=paste(gene, collapse=",")),
                        list(variant.id=id)]
  
  snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
  snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$gene,"\\,")[[1]]]
  
  ### get frequencies
  message("Allele Freqs")
  
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")
  
  af <- data.table(ad=expand.grid(ad$data)[,1],
                   dp=expand.grid(dp)[,1],
                   sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                   variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))
  
  ### tack them together
  message("merge")
  afi <- merge(af, snp.dt1.an, by="variant.id")
  afi <- merge(afi, snp.dt, by.x="variant.id", by.y="id")
  
  afi[,af:=ad/dp]
  
  ### calculate effective read-depth
  afis <- merge(afi, samps, by="sampleId")
  
  afis[chr=="X", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
  afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
  afis[,af_nEff:=round(af*nEff)/nEff]
  
  ### return
  afis[,-c("n"), with=F]
}

### test
data <- getData()

glm.out.padjust %>%
filter(pos %in% haplotags.snps) %>%
    mutate(start=pos,
         end=pos) ->
  haplo_tags_SNPids_and_inv
  
####
haplotag_snps_AFS = foreach(i=1:dim(haplo_tags_SNPids_and_inv)[1], .combine = "rbind", .errorhandling = "remove")%do%{
  
  snp_tmp <- getData(chr="2L", 
                     start=haplo_tags_SNPids_and_inv$pos[i], 
                     end=haplo_tags_SNPids_and_inv$pos[i]) %>%
    mutate(win = haplo_tags_SNPids_and_inv$Win[i]  )
  return(snp_tmp)
}

haplotag_snps_AFS$win %>% table

#### Polarize
haplotag_snps_AFS %>%
  filter(sampleId == "SIM") %>%
  dplyr::select(variant.id, pos, sim_af=af) ->
  SIM_AF

left_join(haplotag_snps_AFS, SIM_AF) %>%
  mutate(af_polarized = case_when(
    sim_af == 0 ~ 1-as.numeric(af),
    sim_af == 1 ~ as.numeric(af))) %>%
  mutate(collectionDate = as.Date(collectionDate, format = "%m/%d/%Y")) ->
  haplotag_snps_AFS_pol

###
haplotag_snps_AFS_pol %>%
  filter(set %in% c("CvilleSet") ) %>%
  filter(!is.na(af_polarized)) %>%
  ggplot(aes(
    x=yday,
    y=af_polarized,
    color = as.factor(variant.id)
  )) +
  geom_line() +
  facet_grid(year~win, scales = "free_x") +
  theme(legend.position = "none") ->
  af_trajectories.pol

ggsave(af_trajectories.pol, file = "af_trajectories.pol.pdf")

###
haplotag_snps_AFS_pol %>%
  filter(set %in% c("CvilleSet") ) %>%
  filter(!is.na(af_polarized)) %>%
  ggplot(aes(
    x=yday,
    y=af,
    color = as.factor(variant.id)
  )) +
  geom_line() +
  facet_grid(year~win, scales = "free_x") +
  theme(legend.position = "none") ->
  af_trajectories.af

ggsave(af_trajectories.af, file = "af_trajectories.af.pdf")

###
save(haplotag_snps_AFS_pol, file = "haplotag_snps_AFS_pol.Rdata")

load("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/haplotag_snps_AFS_pol.Rdata")


######
load("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/nasa_power.weather.mod.Rdata")
names(weather.ave)[1] = "sampleId"

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

weather.ave %>%
  filter(mod == 2) %>%
  dplyr::select(sampleId, temp.max) -> tmax

weather.ave %>%
  filter(mod == 9) %>%
  dplyr::select(sampleId, temp.ave) -> tave

weather.ave %>%
  filter(mod == 8) %>%
  dplyr::select(sampleId, humidity.ave) -> have


cbind(tmax,  tave[,-1], have[,-1]) %>%
  as.data.frame() %>%
  group_by(sampleId) %>%
  slice_head(n=1) -> obs_for_cors


###
haplotag_snps_AFS_pol %>%
  group_by(sampleId, collectionDate, set, year, win, yday) %>%
  filter(set %in% c("CvilleSet", "DrosEU", "DrosRTEC") ) %>%
  filter(!is.na(af_polarized)) %>%
  summarise(Mean_haplotag = mean(af, na.rm = T),
            #sd_haplotag = sd(af, na.rm = T),
            #ci_l = ci(af, na.rm = T)[2],
            #ci_h = ci(af, na.rm = T)[3],
  ) %>%
  left_join(obs_for_cors) %>% 
  filter(set == "CvilleSet") %>%
  filter(!is.na(win)) %>%
  as.data.frame() ->
  Cville_haplotags_for_viz

#####
Cville_haplotags_for_viz %>%
  ggplot(aes(
    x=yday,
    y=Mean_haplotag,
    #ymin=ci_l,
    #ymax=ci_h,
    color=(temp.max),
  )) + 
  #geom_smooth(method = "lm", se = F, size = 0.8, color = "grey") +
  #geom_errorbar(width = 0.1) +
  scale_color_gradient2(low="steelblue", high = "firebrick2", mid = "gold1", 
                        midpoint = 25) +
  geom_point(aes(shape=as.factor(year))) +
  geom_smooth( se = F, size = 0.8, color = "black", linetype = "dashed") +
  ylim(0,0.38) +
  theme_bw() + 
  facet_grid(win~.)->
  haplo.time.colortemp.ave

ggsave(haplo.time.colortemp.ave, file ="haplo.time.colortemp.ave.pdf", h = 6, w = 3.2)



#####
Cville_haplotags_for_viz %>%
  ggplot(aes(
    x=yday,
    y=Mean_haplotag,
    #ymin=ci_l,
    #ymax=ci_h,
    color=(temp.ave),
  )) + 
  #geom_smooth(method = "lm", se = F, size = 0.8, color = "grey") +
  #geom_errorbar(width = 0.1) +
  scale_color_gradient2(low="steelblue", high = "firebrick2", mid = "gold1", 
                        midpoint = 25) +
  geom_point(aes(shape=as.factor(year))) +
  geom_smooth( se = F, size = 0.8, color = "black", linetype = "dashed") +
  ylim(0,0.38) +
  theme_bw() + 
  facet_grid(win~.)->
  haplo.time.color.tempave9

ggsave(haplo.time.color.tempave9, file ="haplo.time.color.tempave9.pdf", h = 6, w = 3.2)

#####
Cville_haplotags_for_viz %>%
  ggplot(aes(
    x=yday,
    y=Mean_haplotag,
    #ymin=ci_l,
    #ymax=ci_h,
    color=(humidity.ave),
  )) + 
  #geom_smooth(method = "lm", se = F, size = 0.8, color = "grey") +
  #geom_errorbar(width = 0.1) +
  scale_color_gradient2(low="steelblue", high = "purple", mid = "grey",
                        midpoint = 75) +
  geom_point(aes(shape=as.factor(year))) +
  geom_smooth( se = F, size = 0.8, color = "black", linetype = "dashed") +
  ylim(0,0.38) +
  theme_bw() + 
  facet_grid(win~.)->
  haplo.time.color.hum8

ggsave(haplo.time.color.hum8, file ="haplo.time.color.hum8.pdf", h = 6, w = 3.2)


#### Plot Ecovariables
Cville_haplotags_for_viz %>% 
  melt(id = c("sampleId", "collectionDate", "set", "year", "win", "yday", "Mean_haplotag")) %>% 
  ggplot(aes(
    x=Mean_haplotag,
    y=value
  )) +
    geom_smooth(
      #method = "lm"
      ) +
    geom_point(aes(    shape=as.factor(year))) +
    facet_wrap(win~variable, scales = "free_y", ncol = 3) ->
    eco_vars_plot
  
  ggsave(eco_vars_plot, file ="eco_vars_plot.pdf", h = 12, w = 8)
  


