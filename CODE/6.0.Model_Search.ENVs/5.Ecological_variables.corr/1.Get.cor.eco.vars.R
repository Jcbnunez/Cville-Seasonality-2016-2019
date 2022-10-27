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
library(broom)


load("/scratch/yey2sn/old_scra/Overwintering_ms/16.Haplotypes/nasa_power.weather.mod.Rdata")
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

cor.test(~ temp.max + temp.ave, dat = obs_for_cors)
cor.test(~ temp.max + humidity.ave, dat = obs_for_cors)
cor.test(~ temp.ave + humidity.ave, dat = obs_for_cors)

#### Plot Top SNPs
#### Plot Top SNPs
#### Plot Top SNPs
#### Plot Top SNPs
#### Plot Top SNPs
#### Plot Top SNPs
#### Plot Top SNPs
#### Plot Top SNPs
#### Plot Top SNPs

load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM_2ltpred_model_and_Files.Rdata")
head(final_in2Lt_markers)


### generate get data
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


####
snps.tracks = foreach(model = c("temp.ave;9;3.Europe_E", "humidity.ave;8;1.Europe_W", "temp.max;2;5.Cville" ), 
        .combine = "rbind")%do%{
  
          message(model)
#### Trajectory analysis
#### base files
base <- "/project/berglandlab/alan/environmental_ombibus_global"
file.tmp <- paste(base, model, paste( model,"glmRNP.Rdata", sep = ".") , sep = "/" )
out.glm <- get(load(file.tmp))

out.glm %>%
  filter(perm == 0) %>%
  filter(chr == "2L") %>% 
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
mutate(win = case_when(
  pos/1e6 > 3.1-0.2 & pos/1e6  < 3.1+0.2 ~ "win_3.1",
  pos/1e6  > 4.7-0.2 & pos/1e6  < 4.7+0.2 ~ "win_4.7",
  pos/1e6  > 5.1-0.2 & pos/1e6  < 5.1+0.2 ~ "win_5.2",
  pos/1e6  > 6.1-0.2 & pos/1e6  < 6.1+0.2 ~ "win_6.1",
  pos/1e6  > 6.8-0.2 & pos/1e6  < 6.8+0.2 ~ "win_6.8",
  pos/1e6  > 9.6-0.2 & pos/1e6  < 9.6+0.2 ~ "win_9.6",
  pos < 3e6 & SNP_id %in% final_in2Lt_markers ~ "left",
  pos > 11e6 & SNP_id %in% final_in2Lt_markers ~ "right")) %>%
  filter(!is.na(win)) %>%
  group_by(win) %>% 
  slice_min(rnp, n = 20) -> snps.of.interest

return(snps.of.interest)
}

snps.tracks %<>%
group_by(win) %>%
  slice_min(rnp, n= 5)
  
#### GET AFS.
#### 
clusts <- fread("./pop.clusters.txt")
snp.traj = 
foreach(snp.i = 1:dim(snps.tracks)[1], .combine = "rbind" )%do%{
  
  sta = snps.tracks$pos[snp.i]
  end = snps.tracks$pos[snp.i]
  mod = snps.tracks$mod[snp.i]
  var = snps.tracks$variable[snp.i]
  
  message(paste(snp.i, sta, mod, var, sep = "|"))
  
  getData(chr="2L", start=sta, end=end) %>%
    mutate(mod = mod, variable = var) %>%
    full_join(clusts)
  
}

snp.traj %>%
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
  mutate(win = case_when(
    pos/1e6 > 3.1-0.2 & pos/1e6  < 3.1+0.2 ~ "win_3.1",
    pos/1e6  > 4.7-0.2 & pos/1e6  < 4.7+0.2 ~ "win_4.7",
    pos/1e6  > 5.1-0.2 & pos/1e6  < 5.1+0.2 ~ "win_5.2",
    pos/1e6  > 6.1-0.2 & pos/1e6  < 6.1+0.2 ~ "win_6.1",
    pos/1e6  > 6.8-0.2 & pos/1e6  < 6.8+0.2 ~ "win_6.8",
    pos/1e6  > 9.6-0.2 & pos/1e6  < 9.6+0.2 ~ "win_9.6",
    pos < 3e6 & SNP_id %in% final_in2Lt_markers ~ "left",
    pos > 11e6 & SNP_id %in% final_in2Lt_markers ~ "right")) %>%
  left_join(obs_for_cors)->
  snp.traj.cands

#### MErge with dtaa


snp.traj.cands %>%
  #filter(year != 2013) %>%
  filter(af > 0.01 & af < 0.99) %>%
  filter(Continental_clusters == "3.Europe_E") %>%
  filter(!is.na(win)) %>%
  ggplot(aes(
    x=temp.ave,
    y= af
  )) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  xlab("Temperature average (45-75 d)") +
  ylab("Allele frequency") + 
  theme_bw() +
  facet_grid(.~win) ->
  plot.EUE

ggsave(plot.EUE, file = "plot.EUE.pdf", w = 8, h = 2.0)


snp.traj.cands %>%
  #filter(year != 2013) %>%
  filter(af > 0.01 & af < 0.99) %>%
  filter(Continental_clusters == "1.Europe_W") %>%
  filter(!is.na(win)) %>%
  ggplot(aes(
    x=humidity.ave,
    y= af
  )) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  xlab("Humidity average (15-45 d)") +
  ylab("Allele frequency") + 
  theme_bw() +
  facet_grid(.~win) ->
  plot.EUW

ggsave(plot.EUW, file = "plot.EUW.pdf", w = 8, h = 2.0)


