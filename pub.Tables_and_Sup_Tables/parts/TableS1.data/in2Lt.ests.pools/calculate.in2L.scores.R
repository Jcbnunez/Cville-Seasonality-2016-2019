### R Script for talk and paper
### 
## This script makes figure 3 of the paper

ingds="/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/PoolSNP_noRep_filter/dest.all.PoolSNP.001.50.10Mar2021.ann.noRep.gds"
inmeta="/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv"
selected_chrs=c("2L")


# Load R packages

library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(tidyverse)
library(gmodels)
library(reshape2)
library(magrittr)

### open GDS file
genofile <- seqOpen(ingds)

### get target populations
samps <- fread(inmeta)

### Include DEST sets
samps <- rbind(
               samps[set=="CvilleSet"]
)

### get subsample of data to work on
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))

## choose number of alleles
snps.dt <- snps.dt[nAlleles==2]



seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt$variant.id)

snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

### select sites
seqSetFilter(genofile, sample.id=samps$sampleId,
             snps.dt[chr%in%selected_chrs]$variant.id)

### get allele frequency data
ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")

dat <- ad$data/dp #if poolSNP
#dat <- ad/dp #if SNAPE
dim(dat)  



## Add metadata
colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")

#Add metadata ad
colnames(ad$data) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(ad$data) <- seqGetData(genofile, "sample.id")

#Add metadata dp
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")

#Generate metadata
left_join(data.frame(sampleId=rownames(dat)), as.data.frame(samps)) -> DEST_DGN_metadata

#### 
dat %>% t() %>% as.data.frame() -> dat.t
dat.t %<>%
  mutate(SNP_id = rownames(.)) %>%
  separate(SNP_id, remove = F, into = c("chr" , "pos"), sep = "_")

dat.t %>%
  melt(id = c("SNP_id", "chr",  "pos")) -> dat.t.m
dat.t.m$pos = as.numeric(dat.t.m$pos)

##Import the inversion mapping SNPs
inversions = read.table(
  "/project/berglandlab/Dmel_genomic_resources/Inversions/inversion_makers_kapun.txt",
  head = T)

names(inversions)[2:3] = c("chr","pos")

###
left_join(inversions, 
          dat.t.m, 
          by = c("chr","pos")) %>%
  .[complete.cases(.),] %>% 
  group_by(inversion,variable) %>%
  summarize(Inv_freq = mean(value)) %>% 
  dcast(variable~inversion, value.var = "Inv_freq") ->
  Inversion_frequencies

write.table(Inversion_frequencies, file = "pools.Inversion_frequencies.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")



####


rm(list = ls())

#load data
# This R object was premade in script 1 ==> 1.Import_GDStoR.r
data_in="/scratch/yey2sn/old_scra/Overwintering_ms/1.Make_Robjects_Analysis/DEST.2.0.poolSNP.Spatial.Temporal.ECfiltered.Rdata"

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

##Import the inversion mapping SNPs
inversions = read.table(
  "/project/berglandlab/Dmel_genomic_resources/Inversions/inversion_makers_kapun.txt",
  head = T)

names(inversions)[2:3] = c("chr","pos")

## Load the genotype matrix
load(data_in)

#### Begin Measure inversions
### Measure inversion proportions
inversions %>% head

dat_filtered_t %>% 
  separate(SNP_id, 
           into = c("chr","pos"), 
           sep = "_", 
           remove = F) ->
  dat_filtered_t_chr_pos

dat_filtered_t_chr_pos$pos = as.numeric(dat_filtered_t_chr_pos$pos)

left_join(inversions, 
          dat_filtered_t_chr_pos, 
          by = c("chr","pos")) %>%
  .[complete.cases(.),] %>%
  melt(id = c("inversion", 
              "chr", 
              "pos", 
              "SNP_id", 
              "allele")) %>% 
  group_by(inversion,variable) %>%
  summarize(Inv_freq = mean(value)) %>% 
  dcast(variable~inversion, value.var = "Inv_freq") ->
  Inversion_frequencies
