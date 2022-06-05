rm(list = ls())

### libraries
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(magrittr)
#use doParallel package to register multiple cores that can be used to run loops in parallel
registerDoParallel(4)

### feed in user arguments
args = commandArgs(trailingOnly=TRUE)
g = as.numeric(args[1])-1
### Launch the array 1-101 <<----- 

###
out_file = "/scratch/yey2sn/Overwintering_ms/4.2.env.omibus.mods/out_folder/"

##load SNP file
load("/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt_25percMissing.Rdata")

#####
##ingds="/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds"
##inmeta="/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv"
###outfile="DEST.2.0.poolSNP.Spatial.Temporal.AllDat"
##
##samps <- fread(inmeta)
##
##### open GDS file
##genofile <- seqOpen(ingds)
##
##### get subsample of data to work on
##seqResetFilter(genofile)
##seqSetFilter(genofile, sample.id=samps$sampleId)
##
##snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
##                      pos=seqGetData(genofile, "position"),
##                      variant.id=seqGetData(genofile, "variant.id"),
##                      nAlleles=seqNumAllele(genofile),
##                      missing=seqMissing(genofile, .progress=T))
##
##snps.dt %>% group_by(chr) %>% summarize(N=n())

#### Generate sets of time models
#sets <- data.table(mod=c(-1, 0, 1:11),
#                   label=LETTERS[1:13],
#                   year=c(NA, rep(1, 12)),
#                   start=c(NA, NA, 0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
#                   end=	 c(NA, NA, 7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))


#library(ggplot2)
#make a list of paths to Alan's files
fl <- paste0("/project/berglandlab/alan/environmental_ombibus/rawData/job", c(1:5000), ".Rdata")
#use a for loop to save file for each permutation

#foreach(g = c(0:100))%do% { 

out = foreach(f = c(1:5000), 
              .combine = "rbind") %do% {
                
                #f = 2
                message(f)
                file = fl[f]
                print(file)
                load(file)
                #sort out observed. 
                glm.out = glm.out[perm == g]
                glm.out
                
                names(glm.out)[9] = "id"
                
                #names(glm.out) = c("mod","AIC","b_temp","se_temp","nObs","p_lrt","variant.id","perm")
                
                left_join(glm.out, sets) %>%
                  left_join(snp.dt[,c("id","chr", "pos")]) %>% 
                  mutate(job = f) ->
                  glm.out.annot
                
                return(glm.out.annot)
                
              }

save(out, 
     file = paste(out_file, 
                  "VA_mod.newmods." , g, ".Rdata", sep = "")
) ## close save

#}
