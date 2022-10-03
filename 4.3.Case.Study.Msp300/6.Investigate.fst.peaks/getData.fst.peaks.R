#### find dirver mutations in w5.2
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
library(broom)

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
#
#CHROM BIN_START BIN_END N_VARIANTS WEIGHTED_FST MEAN_FST
#2L      5170001 5180000         20        0.803    0.537
#2L      5190001 5200000        200        0.719    0.595
#2L      5185001 5195000        212        0.706    0.565
#2L      5200001 5210000        205        0.690    0.546

data.win1 <- getData(chr="2L", start=5170001, end=5180000) %>%
  mutate(peak="5170001_5180000") %>%
  group_by(variant.id) %>%
  dplyr::select(chr,pos,col, gene) %>%
  slice_head()

data.win2 <- getData(chr="2L", start=5190001, end=5200000) %>%
  mutate(peak="5190001_5200000") %>%
  group_by(variant.id) %>%
  dplyr::select(chr,pos,col, gene) %>%
  slice_head()

data.win3 <- getData(chr="2L", start=5185001, end=5195000) %>%
  mutate(peak="5185001_5195000") %>%
  group_by(variant.id) %>%
  dplyr::select(chr,pos,col, gene) %>%
  slice_head()

data.win4 <- getData(chr="2L", start=5200001, end=5210000) %>%
  mutate(peak="5200001_5210000") %>%
  group_by(variant.id) %>%
  dplyr::select(chr,pos,col, gene) %>%
  slice_head()

rbind(data.win1, data.win2, data.win3, data.win4) ->
  win.snp.annots

###
glm.cville <- get(load("/project/berglandlab/alan/environmental_ombibus_global/temp.max;2;5.Cville/temp.max;2;5.Cville.glmRNP.Rdata"))
glm.cville.0 = glm.cville %>%
  filter(perm == 0)

left_join(glm.cville.0, win.snp.annots) %>%
  filter(col == "missense_variant") %>%
  filter(rnp < 0.05) %>%
  group_by(pos) %>%
  slice_head() %>% 
  mutate(tsp.mel = 5192177,
         tsp.sim = 4996892) %>%
  as.data.frame() %>%
  mutate(diff.to.tsp = 5192177-pos) %>%
  mutate(sim.pos = tsp.sim- diff.to.tsp)

