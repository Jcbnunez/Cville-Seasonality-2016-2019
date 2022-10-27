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
###
###### Summarize trajectories over the haplotypes
setwd("/scratch/yey2sn/Overwintering_ms/12.trajectory_analysis/")
## We will use generated in step 4
load("./haplo_tags_SNPids.Rdata")
haplo_tags_SNPids

  haplo_tags_SNPids %<>% 
  filter(class == "GLM_LD") %>%
  separate(SNP_id, into = c("chr","pos", "type"), remove = F) %>%
  mutate(win = 
           case_when(
             pos > 4650065 & pos < 4799922 ~ "win_4.7",
             pos > 5100324 & pos < 5349218 ~ "win_5.2",
             pos > 6100321 & pos < 6349489 ~ "win_6.2",
             pos > 9500286 & pos < 9700005 ~ "win_9.6"
           ))  

load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM_2ltpred_model_and_Files.Rdata")
final_in2Lt_markers %<>%
  data.frame(SNP_id = .) %>%
  separate(SNP_id, into = c("chr", "pos", "type"), sep = "_", remove = F) %>%
  mutate(win = "inv",
         class = "inv")

rbind(haplo_tags_SNPids, final_in2Lt_markers) -> haplo_tags_SNPids_and_inv

#save(haplo_tags_SNPids_and_inv , file = "/project/berglandlab/jcbnunez/Shared_w_Alan/haplo_tags_SNPids.Rdata")

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


###
### GET AFS
###
haplotag_snps_AFS = foreach(i=1:dim(haplo_tags_SNPids_and_inv)[1], .combine = "rbind")%do%{
  
  snp_tmp <- getData(chr="2L", 
                     start=haplo_tags_SNPids_and_inv$pos[i], 
                     end=haplo_tags_SNPids_and_inv$pos[i]) %>%
    mutate(win = haplo_tags_SNPids_and_inv$win[i]  )
    return(snp_tmp)
}

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


### add temp
load("/scratch/yey2sn/Overwintering_ms/4.2.env.omibus.mods/weather.Rdata")
names(weather.ave)[1] = "sampleId"

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

weather.ave %>%
  filter(mod == 2) %>%
  dplyr::select(sampleId, temp.max) -> tmax

weather.ave %>%
  filter(mod == 3) %>%
  dplyr::select(sampleId, precip.var) -> pvar

weather.ave %>%
  filter(mod == 4) %>%
  dplyr::select(sampleId, temp.ave) -> tave

weather.ave %>%
  filter(mod == 11) %>%
  dplyr::select(sampleId, precip.ave) -> pave

cbind(tmax, pvar[,-1], tave[,-1], pave[,-1]) %>%
  as.data.frame() %>%
  group_by(sampleId) %>%
  slice_head(n=1) -> obs_for_cors

#weather.ave %>%
#  left_join(sets) %>%
#  dplyr::select(temp.ave, temp.max, precip.ave)
  
###
save(haplotag_snps_AFS_pol, file = "haplotag_snps_AFS_pol.Rdata")

load("haplotag_snps_AFS_pol.Rdata")
###
haplotag_snps_AFS_pol %>%
  filter(set %in% c("CvilleSet", "DrosEU", "DrosRTEC") ) %>%
  filter(!is.na(af_polarized)) %>%
  group_by(sampleId, collectionDate, set, year, win, yday) %>%
  summarise(Mean_haplotag = mean(af, na.rm = T),
            #sd_haplotag = sd(af, na.rm = T),
            #ci_l = ci(af, na.rm = T)[2],
            #ci_h = ci(af, na.rm = T)[3],
            ) %>%
  filter(set == "CvilleSet") %>% 
  left_join(obs_for_cors) %>% 
  as.data.frame() ->
  Cville_haplotags_for_viz
  
#####
Cville_haplotags_for_viz %>%
  ggplot(aes(
    x=yday,
    y=Mean_haplotag,
    ymin=ci_l,
    ymax=ci_h,
    color=fahrenheit.to.celsius(temp.max),
  )) + 
  #geom_smooth(method = "lm", se = F, size = 0.8, color = "grey") +
  geom_errorbar(width = 0.1) +
  scale_color_gradient2(low="steelblue", high = "firebrick2", mid = "gold1", 
                        midpoint = 25) +
  geom_point(aes(shape=as.factor(year))) +
  geom_smooth( se = F, size = 0.8, color = "black", linetype = "dashed") +
  ylim(0,0.38) +
  theme_bw() + 
  facet_grid(win~.)->
  haplo.time.colortemp.ave
ggsave(haplo.time.colortemp.ave, file ="haplo.time.colortemp.ave.pdf", h = 6, w = 4.5)



Cville_haplotags_for_viz %>%
  ggplot(aes(
    x=precip.ave,
    y=Mean_haplotag,
    ymin=ci_l,
    ymax=ci_h,
    color=as.factor(year),
  )) + 
  geom_errorbar(width = 0.001) +
  #scale_color_gradient2(low="blue", high = "red", mid = "yellow", midpoint = 13) +
  geom_point(aes(shape=as.factor(year))) +
  geom_smooth( se = F, size = 0.8, #color = "black", 
               linetype = "dashed", method = "lm") +
  theme_bw() + 
  ylim(0,0.38) +
  theme(legend.position = "none") +
  ggtitle("Precipt Avg") +
  facet_grid(win~year)->
  haplo_mean_pcip
ggsave(haplo_mean_pcip, file ="haplo_mean_pcip.pdf", h = 6, w = 5)


Cville_haplotags_for_viz %>%
  ggplot(aes(
    x=precip.var,
    y=Mean_haplotag,
    ymin=ci_l,
    ymax=ci_h,
    color=as.factor(year),
  )) + 
  geom_errorbar(width = 0.001) +
  #scale_color_gradient2(low="blue", high = "red", mid = "yellow", midpoint = 13) +
  geom_point(aes(shape=as.factor(year))) +
  geom_smooth( se = F, size = 0.8, #color = "black", 
               linetype = "dashed", method = "lm") +
  theme_bw() + 
  ylim(0,0.38) +
  theme(legend.position = "none") +
  ggtitle("precip.var") +
  facet_grid(win~year)->
  precip.var_plot
ggsave(precip.var_plot, file ="precip.var_plot.pdf", h = 6, w = 5)


Cville_haplotags_for_viz %>%
  ggplot(aes(
    x=temp.ave,
    y=Mean_haplotag,
    ymin=ci_l,
    ymax=ci_h,
    color=as.factor(year),
  )) + 
  geom_errorbar(width = 0.1) +
  #scale_color_gradient2(low="blue", high = "red", mid = "yellow", midpoint = 13) +
  geom_point(aes(shape=as.factor(year))) +
  geom_smooth( se = F, size = 0.8, #color = "black", 
               linetype = "dashed", method = "lm") +
  theme_bw() + 
  ylim(0,0.38) +
  theme(legend.position = "none") +
  ggtitle("T Ave") +
facet_grid(win~year) ->
  haplo_T_ave
ggsave(haplo_T_ave, file ="haplo_T_ave.pdf", h = 6, w = 5)




ggsave(haplo_mean_temp + haplo_mean_time, file ="Cville_haplotag_trajectories.pdf", w = 5.5, h =6)

