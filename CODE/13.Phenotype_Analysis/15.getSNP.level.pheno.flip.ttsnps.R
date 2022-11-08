###find tt snps dosage
#this script records the genotype of each of our TT loci for each of our DGRP lines
#the original does like 4 other things, that we need to clean up

###find tt snps within each peak - and the phenotypes they correspond to

library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(matrixStats)
library(SeqArray)
library(magrittr)
registerDoParallel(4)
#load in sliding window data

#load in peak and window data
#setwd("/scratch/bal7cg/Deficiency-Line-confirmation/")

dt.coenrich = get(load("/project/berglandlab/jcbnunez/Shared_w_Alan/Junemodels.slidingwindow.Rdata"))

realglm = dt.coenrich %>% 
  select(c(snp.positions, chr, gwas.pheno, start, end, win.i, perm.id)) %>% 
  filter(perm.id== 0)



perm.o.chr = foreach(chr.i = unique(realglm$chr),
                     .combine = "rbind") %do%{ ### open chr.i
                       
                       #message(chr.i)
                       
                       perm.o = foreach(p = unique(realglm$gwas.pheno), 
                                        .combine = "rbind", 
                                        .errorhandling = "remove") %dopar% {   ### opens p
                                          
                                          #p  = unique(realglm$gwas.pheno)[17]
                                          #chr.i = "3R"
                                          #
                                          message(paste(p, chr.i, sep = "|"))
                                          
                                          zoom = realglm[gwas.pheno == p][chr == chr.i]
                                          
                                          #different peakstart/ end if trying to find ttsnps for a certain gene
                                          #range of msp300 is 5100877-5207002
                                          # peakstart = 5100877
                                          # peaksend = 5207002
                                          
                                          # if(dim(zoom)[1] > 0){
                                          zoom = zoom %>% 
                                          select(snp.positions, gwas.pheno)
                                          
                                          #zoom$pasted = paste0(zoom$group, zoom$snp.positions)
                                          split = strsplit(zoom$snp.positions, split = ",")
                                          table = plyr::ldply(split, rbind)
                                          full = cbind(zoom, table)
                                          full = full[,-1]
                                          #what are our phenotnames?
                                          
                                          melt = melt(full, id.vars = c("gwas.pheno"))
                                          #remove meaninglyless "variable" column
                                          melt = melt[, -2]
                                          names(melt)[which(names(melt) == "value")] = "position"
                                          melt[,chr:=chr.i]
                                          melt = melt[complete.cases(melt$position)]
                                          return(melt) 
                                        } ### closes p 
                       
                       return(perm.o)
                       
                     } #### closes chri.i

perm.o.chr %<>%
  mutate(SNP_id = paste(chr, position, sep = "_"))

#perm.bind = rbindlist(perm.o)
#now turn back around, and create a datable that lists phenotypes as a vector for each tt snp
ttsnps = unique(perm.o.chr$SNP_id)

tt.out = foreach(g = ttsnps,
                 .combine = "rbind"
                 )%do% { ## open g
  message(g)
 # g = ttsnps[1]
  ttdata = perm.o.chr[SNP_id == g]
  phenos = ttdata$gwas.pheno
  chr = unique(ttdata$chr)
  pos = unique(ttdata$position)
  
  #create a row of what will ultimately our snp data table.
  o = data.frame ( 
    chr = chr,
    pos = pos,
    SNP_id = g,
    N.phenos = length(phenos),
    Description = toString(phenos)
  )
  
  return(o)
                 } ### close g


#ttbind$GWAS_sig = "TT"
save(tt.out, file = "SNP.phenos.Rdata")
#write.csv(ttbind, "/home/bal7cg/SNP.table.txt")
