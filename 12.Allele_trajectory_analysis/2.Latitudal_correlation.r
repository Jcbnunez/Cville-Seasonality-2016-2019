### libraries

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
library(vcfR)
library(scatterpie)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
registerDoMC(2)

####
setwd("/scratch/yey2sn/Overwintering_ms/12.trajectory_analysis/")

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

#######
#######
#######
#######
glm.file <- "/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_0.Rdata"
load(glm.file)
glm.out %>%
  filter(mod == "aveTemp+year_factor",
         chr == "2L",
         rnp.clean < 0.05) -> outliers

load("/project/berglandlab/DEST_Charlottesville_TYS/weatherAve.Rdata")
names(weather.ave)[1] = "sampleId"

padding=0
starts= c(5155762, 6255762, 9505762)
ends= c(5255762, 6355762, 9605762)

win_5np = getData(chr="2L", start=5155762-padding , end=5255762+padding) 
win_6np = getData(chr="2L", start=6255762-padding , end=6355762+padding) 
win_9np = getData(chr="2L", start=9505762-padding , end=9605762+padding)

### join datasets
rbind(mutate(win_5np, win = "1.win5"), 
      mutate(win_6np, win = "2.win6"),
      mutate(win_9np, win = "3.win9")) %>%   
  left_join(weather.ave) ->
  win_outliers_no_padding


### Make additional files;, simulans polarization and annotations
win_outliers_no_padding %>%
  filter( set %in% c("dgn"),
          sampleId == "SIM") %>%
  dplyr::select(chr, pos, sim_af = af) -> sim_af

win_outliers_no_padding %>%
  group_by(pos) %>%
  slice_head(n=1) %>%
  dplyr::select(pos, col) -> annots

#######
#######
#######
#######

win_outliers_no_padding %>%
  left_join(sim_af) %>%
  #filter(pos %in% outliers$pos ) %>%
  mutate(af_neff_pol = case_when(
    sim_af == 0 ~ 1-af_nEff,
    sim_af == 1 ~ af_nEff))  -> 
  win_outliers_no_padding_DEST_polarized


### list of targets
win_outliers_no_padding_DEST_polarized$pos %>% unique -> pos_int
length(pos_int)
### run model
latitude_model =foreach(i=1:length(pos_int), .errorhandling = "remove" , .combine = "rbind")%do%{
  message(round(i/length(pos_int), 3))  
  filter(win_outliers_no_padding_DEST_polarized, 
         pos == pos_int[i], 
         set %in% c("DrosRTEC","DrosEU"),
         season == "spring") -> tmp
  tmp %<>% filter(!is.na(af_neff_pol))
  summary(lm(af_nEff ~ lat, data = tmp)) -> model_tmp
  data.frame(
    pos=pos_int[i],
    lat_beta = model_tmp$coefficients[2,1],
    lat_p_val = model_tmp$coefficients[2,4]
  ) -> model_out
  return(model_out)
}

save(latitude_model, file = "latitude_model.Rdata")
