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

in2Lmarker <- vroom("./inv2L_correlated_markers_Dm3.txt")

in2Lmarker %>%
  filter(cor_rank >= 0.6) ->
  in2Lmarker.0.6


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
in2Lmarker.0.6 %>% 
  dplyr::select(SNP_id, correlation, p.value.bonfe) %>% 
  left_join(snp.dt) ->
  in2Lmarker.0.6.annot

#### make pools of controls
snp.dt %>%
  filter(chr != "2L") %>%
  filter(invName != "none") %>%
  left_join(snp.dt[,c("chr", "pos", "cm_mb", "af")]) ->
  pool_of_controls

#### choose contols

in.df = in2Lmarker.0.6.annot

in.df = in2Lmarker.0.6.annot
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
                              slice_head() %>% 
                              mutate(control_snp = paste(chr, pos, "SNP", sep = "_" )) %>%
                              mutate(matched_to = tmp.achor$SNP_id, type = "control")->
                              select.control
                            
                    return(select.control)
                            
                          }


write.table(chosen_controls, 
            file = "Matched.controls.In2Lt.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

