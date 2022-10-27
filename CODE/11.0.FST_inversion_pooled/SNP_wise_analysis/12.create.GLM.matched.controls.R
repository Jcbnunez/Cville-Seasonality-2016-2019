##### ----> extract the In2Lt markers
##### 
##### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(foreach)
library(doParallel)

####
####

glm.snps <- get(load("./Cville.GLM.rnp5.snps.Rdata"))

### load master snp object
### Load master SNP file
load("/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt_25percMissing.Rdata")
head(snp.dt)
snp.dt %<>%
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_"))


genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)
### get subsample of data to work on
seqResetFilter(genofile)
#seqSetFilter(genofile, sample.id=samps.cville$sampleId)
seqSetFilter(genofile, variant.id=snp.dt$id)
snp.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]


#### make in2Lt markers
#### 
glm.snps %>% 
  left_join(snp.dt) ->
  glm.snps.annot

#### make pools of controls
snp.dt %>%
  filter(chr != "2L") %>%
  filter(invName != "none") %>%
  left_join(snp.dt[,c("chr", "pos", "cm_mb", "af")]) ->
  pool_of_controls

#### choose contols

in.df = glm.snps.annot
cont.df = pool_of_controls

chosen_controls = foreach(i = 1:dim(in.df)[1],
                          .combine = "rbind",
                          .errorhandling ="remove")%do%{
                            
                            message(paste(i, dim(in.df)[1], sep = " | "))
                            in.df[i,] ->
                              tmp.achor
                            
                            cont.df %>%
                              filter(cm_mb > tmp.achor$cm_mb-0.20 & cm_mb < tmp.achor$cm_mb+0.20,
                                     af > tmp.achor$af-0.030 & af < tmp.achor$af+0.030
                              ) %>%
                              slice_sample(n=100) %>% 
                              mutate(control_snp = paste(chr, pos, "SNP", sep = "_" )) %>%
                              mutate(matched_to = tmp.achor$SNP_id, type = "control")->
                              select.control
                            
                            return(select.control)
                            
                          }


write.table(chosen_controls, 
            file = "Matched.controls.quantiles.GLM.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

