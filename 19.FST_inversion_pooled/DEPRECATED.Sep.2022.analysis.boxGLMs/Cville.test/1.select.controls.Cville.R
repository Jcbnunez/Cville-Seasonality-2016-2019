#### FST for the inversion regions
#### FST for the inversion regions
#### 
#### 
#### Script 1. Select the control markers 
#### 
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

### Load master SNP file
load("/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt_25percMissing.Rdata")

head(snp.dt)

#### load model output
base <- "/project/berglandlab/alan/environmental_ombibus_global"
models = c("temp.max;2;5.Cville")
k=1
file <- paste(base, models[k], paste(models[k],"glmRNP.Rdata", sep = ".") , sep = "/" )
print(file)

message(models[k])
out.glm <- get(load(file))

### filter
out.glm %>%
  filter(perm == 0) %>%
  filter(rnp < 0.05) %>%
  filter(chr == "2L") %>% 
  filter(invName == "2Lt") ->
  in2lt.outliers

#### import genomic data
genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)

samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
samps[nchar(month)==1, month:=paste("0", month, sep="")]
samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
samps[nchar(day)==1, day:=paste("0", day, sep="")]
samps <- samps[set!="dgn"]
samps[,Date:=date(paste(year, month, day, sep="-"))]

###
samps %>%
  filter(set == "CvilleSet") ->
  samps.cville

### get subsample of data to work on
seqResetFilter(genofile)
#seqSetFilter(genofile, sample.id=samps.cville$sampleId)
seqSetFilter(genofile, sample.id=samps.cville$sampleId, variant.id=snp.dt$id)
snp.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

#### generate matched controls
snp.dt %>%
  filter(id %in% in2lt.outliers$variant.id) -> glm.snps

snp.dt %>%
  filter(chr != "2L" & chr != "X") %>%
  filter(VA_ch == TRUE) %>%
  filter(is.na(libs))-> pool.of.controls

##### Samplinc scheme
macthed.controls = 
  foreach(i = 1:dim(glm.snps)[1], 
        .combine = "rbind",
        .errorhandling = "remove")%do%{ ## open i
          
          message(i)
          
          glm.snps[i,] -> target.tmp
          
          pool.of.controls %>%
            filter(af >= target.tmp$af - 0.01 & af <= target.tmp$af + 0.01  ) %>%
            filter(cm_mb >= target.tmp$cm_mb - 0. & cm_mb <= target.tmp$cm_mb + 0.1   ) %>%
            filter(invName == "none") %>%
            sample_n(1) -> matched.noInv
          
          pool.of.controls %>%
            filter(af >= target.tmp$af - 0.01 & af <= target.tmp$af + 0.01  ) %>%
            filter(cm_mb >= target.tmp$cm_mb - 0. & cm_mb <= target.tmp$cm_mb + 0.1   ) %>%
            filter(invName != "none") %>%
            sample_n(1) -> matched.Inv
          
          rbind(mutate(matched.noInv, matched.type = "noInv", 
                       matched.to = paste(target.tmp$id,target.tmp$chr,target.tmp$pos, sep = "_" ) ),
                     mutate(matched.Inv, matched.type = "Inv", 
                            matched.to = paste(target.tmp$id,target.tmp$chr,target.tmp$pos, sep = "_" ))
                     )
          
        } ### close i

macthed.controls %>% 
  filter(matched.type == "noInv") %>%
  group_by(id) %>%
  slice_head() ->
  macthed.controls.noInv
  
macthed.controls %>% 
  filter(matched.type == "Inv") %>%
  group_by(id) %>%
  slice_head() ->
  macthed.controls.Inv

########### save all
save(
  glm.snps,
  macthed.controls.noInv,
  macthed.controls.Inv,
  file = "glm.and.matchedControls.Rdata"
)


