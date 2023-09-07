
#args = commandArgs(trailingOnly=TRUE)
#job=as.numeric(args[1])-1
#k=as.numeric(args[1])-1
#model=args[2]
#k=0


### libraries
library(data.table)
library(tidyverse)
library(magrittr)
library(foreach)
library(doMC)
registerDoMC(4)
library(tidyr)

####

####  create window objects
final.windows.pos = 
  data.frame(win.name = c("left", "win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6", "right" ),
             mid = c(2.2, 3.0, 4.67, 5.12, 6.2, 6.8 , 9.6, 13.1),
             chr = "2L"
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )

### load thermal GLM object
sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

models = c(
		   #"temp.max;4;5.Cville",
		   "humidity.ave;10;1.Europe_W",
		   "humidity.var;3;3.Europe_E",
		   "temp.min;8;2.North_America_E"
		    #,
           #"temp.ave;9;3.Europe_E" #,
           #"temp.ave;1;2.North_America_E" #,
           #"humidity.ave;8;1.Europe_W"
           )

base <- "/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper"


######
######
### ----> loac Cville object ... this is a constant througt
#file.cvile <- paste(base, "temp.max;2;5.Cville", "temp.max;2;5.Cville.glmRNP.Rdata", sep = "/" )
#print(file.cvile)
file.cvile <- paste(base, "Revision_Best_Models", paste( "temp.max;2;5.Cville","v2.glmRNP.Rdata", sep = ".") , sep = "/" )
  print(file.cvile)

out.glm.cvile <- get(load(file.cvile))



#### load and plot
load("./window.enrich.set.Rdata")

#####
##### ---> Cline analysis
#####
##### ---> Cline analysis
#####
##### ---> Cline analysis

head(out.glm.cvile)


### load core20 data
core20.orig <- fread("/project/berglandlab/alan/drosRTEC/mnt/pricey_1/dropPop/mel_clinal_uniquepops.glm")
#core20.swap <- fread("/project/berglandlab/alan/drosRTEC/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm")

core20.orig[,set:="cline"]
#core20.swap[,set:="swap"]
core20 <- rbind(core20.orig, core20.swap)
core20 <- core20.orig

setnames(core20, "chrom", "chr")

### R5 -> R6 DGRP conversion table
liftover.fn <- "/project/berglandlab/Dmel_genomic_resources/liftOver_files/dest.all.PoolSNP.001.50.dm3.dm6.csv"
liftover <- fread(liftover.fn)
liftover[,SNP:=paste(dm3_chr, dm3_pos, "SNP", sep="_")]

### do liftover
setnames(core20, c("chr", "pos"), c("dm3_chr", "dm3_pos"))
setkey(core20, dm3_chr, dm3_pos)
setkey(liftover, dm3_chr, dm3_pos)

core20.or <- merge(core20, liftover)

setnames(core20.or, c("dm6_chr", "dm6_pos"), c("chr", "pos"))

### merge with VA_glm
#setkey(core20, chr, pos)
#setkey(glm.out, chr, pos)
#m <- merge(glm.out, core20)

    ####
    #out.glm.matched %>%
    #  filter(perm == k) -> out.glm.matched.mod
    out.glm.matched.mod = core20.or
    setnames(out.glm.matched.mod, c("clinal.coef", "clinal.p" , "SNP"), c("b_temp", "p_lrt" ,"snp.id"))
    
    #b_temp
    #se_temp
    
    out.glm.cvile %>%
      filter(perm == k) %>%
      mutate(snp.id = paste(chr, pos, "SNP", sep = "_")) -> out.glm.cvile.mod
    
    ### merge
    message("Merge datasets")
    
    m1 <- merge(out.glm.cvile.mod, 
                out.glm.matched.mod, 
                by.x="snp.id", by.y="snp.id")
    
    m1 <- m1[chr.x!="X"]
    m1[!is.na(p_lrt.x) & !is.na(p_lrt.y) ,cville.rank := rank(p_lrt.x)/length(p_lrt.x)]
    m1[!is.na(p_lrt.x) & !is.na(p_lrt.y) ,anchor.rank := rank(p_lrt.y)/length(p_lrt.y)]
    
    message("Estimating the betas")
    
    m1.beta <- m1[,list(cville.beta=b_temp.x,
                        anchor.beta=b_temp.y), 
                  list(variant.id)]
    
    m1 <- merge(m1, m1.beta, by="variant.id")
    m1[,inv:=invName!="none"]
    
    thrs <- c(
      0.05
    )
    
    wins=final.windows.pos
    model = "cline.Machado"
    
    enrichment.cline <- foreach(#chr.i=unique(m1$chr.x), 
      chr.i="2L",
      .combine="rbind", 
      .errorhandling="remove")%do%{
        
        foreach(pos.ith=1:dim(filter(wins, chr == chr.i))[1], 
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  foreach(thr.i=thrs, 
                          .combine="rbind", 
                          .errorhandling="remove")%do%{
                            
                            # chr.i <- "3R"; inv.i <- T; thr.i<-0.05
                            message(paste(chr.i, pos.ith, dim(filter(wins, chr == chr.i))[1] , thr.i, sep=" / "))
                            
                            ### Machado set, enrichment
                            tab <- table(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ]$cville.rank  < thr.i,
                                         m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ]$anchor.rank <  thr.i)
                            
                            fet <- fisher.test(tab)
                            
                            ### Machado set, sign test
                            st.T <- sum(sign(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ][cville.rank<thr.i & anchor.rank<thr.i]$cville.beta) ==
                                          sign(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ][cville.rank<thr.i & anchor.rank<thr.i]$anchor.beta))
                            st.F <- sum(sign(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ][cville.rank<thr.i & anchor.rank<thr.i]$cville.beta) !=
                                          sign(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ][cville.rank<thr.i & anchor.rank<thr.i]$anchor.beta))
                            
                            bt <- binom.test(st.T, st.T+st.F, .5)
                            
                            
                            tmp1 <- data.table(chr.x=chr.i, 
                                               thr=thr.i , 
                                               perm=k,
                                               anchor.model=model,
                                               win.start=wins$start[pos.ith],
                                               win.end=wins$end[pos.ith],
                                               or=fet$estimate, p=fet$p.value, lci=fet$conf.int[1], uci=fet$conf.int[2],
                                               st.pr=bt$estimate, st.p=bt$p.value, st.lci=bt$conf.int[1], st.uci=bt$conf.int[2])
                            
                            ### return
                            return(tmp1)
                            
                          }
                }
      }
    
    #o.win %>%
    #  filter(p < 0.05)
    
    #####
    ##### ---> Core 20
    #####
    ##### ---> Core 20
    #####
    ##### ---> Core 20 <<< stopped here
    head(out.glm.cvile)
    
    ### load core20 data
    ### load core20 GLM object
    load(file="/project/berglandlab/alan/core20glm.Rdata")
    core20.glm %>%
      filter(mod == "season+locality_factor" & set == "no_va") ->
      core20.orig
  
    #core20.orig <- fread("/project/berglandlab/alan/drosRTEC/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear.f_s.glm")
    #core20.swap <- fread("/project/berglandlab/alan/drosRTEC/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm")
    
    core20.orig[,set:="orig"]
    #core20.swap[,set:="swap"]
    #core20 <- rbind(core20.orig, core20.swap)
    
    core20 <- core20.orig
    #setnames(core20, "chrom", "chr")
    
    #### R5 -> R6 DGRP conversion table
    #liftover.fn <- "/project/berglandlab/Dmel_genomic_resources/liftOver_files/dest.all.PoolSNP.001.50.dm3.dm6.csv"
    #liftover <- fread(liftover.fn)
    #liftover[,SNP:=paste(dm3_chr, dm3_pos, "SNP", sep="_")]
    #
    #### do liftover
    #setnames(core20, c("chr", "pos"), c("dm3_chr", "dm3_pos"))
    #setkey(core20, dm3_chr, dm3_pos)
    #setkey(liftover, dm3_chr, dm3_pos)
    #
    #core20 <- merge(core20, liftover)
    #
    #setnames(core20, c("dm6_chr", "dm6_pos"), c("chr", "pos"))
    
    m1.beta <- core20[,list(seas.beta=tstrsplit(b, ";")%>%last%>%as.numeric), list(variant.id)]
    out.glm.matched.mod = merge(core20, m1.beta)
    out.glm.matched.mod %<>%
      dplyr::select(variant.id, perm, p.lrt, seas.beta)
    
    setnames(out.glm.matched.mod, c("seas.beta", "p.lrt" ), c("b_temp", "p_lrt"))
    
##

    ### merge
    message("Merge datasets")
    
    m1 <- merge(out.glm.cvile.mod, 
                out.glm.matched.mod, 
                by.x="variant.id", by.y="variant.id")
    
    m1 <- m1[chr !="X"]
    m1[!is.na(p_lrt.x) & !is.na(p_lrt.y) ,cville.rank := rank(p_lrt.x)/length(p_lrt.x)]
    m1[!is.na(p_lrt.x) & !is.na(p_lrt.y) ,anchor.rank := rank(p_lrt.y)/length(p_lrt.y)]
    
    message("Estimating the betas")
    
    m1.beta <- m1[,list(cville.beta=b_temp.x,
                        anchor.beta=b_temp.y), 
                  list(variant.id)]
    
    m1 <- merge(m1, m1.beta, by="variant.id")
    m1[,inv:=invName!="none"]
    m1 %<>%
      mutate(chr.x = chr, pos.x = pos)
    
    thrs <- c(
      0.05
    )
    
    wins=final.windows.pos
    model = "core20.sansVA.Machado"
    
    enrichment.core20 <- foreach(#chr.i=unique(m1$chr.x), 
      chr.i="2L",
      .combine="rbind", 
      .errorhandling="remove")%do%{
        
        foreach(pos.ith=1:dim(filter(wins, chr == chr.i))[1], 
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  foreach(thr.i=thrs, 
                          .combine="rbind", 
                          .errorhandling="remove")%do%{
                            
                            # chr.i <- "3R"; inv.i <- T; thr.i<-0.05
                            message(paste(chr.i, pos.ith, dim(filter(wins, chr == chr.i))[1] , thr.i, sep=" / "))
                            
                            ### Machado set, enrichment
                            tab <- table(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ]$cville.rank  < thr.i,
                                         m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ]$anchor.rank <  thr.i)
                            
                            fet <- fisher.test(tab)
                            
                            ### Machado set, sign test
                            st.T <- sum(sign(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ][cville.rank<thr.i & anchor.rank<thr.i]$cville.beta) ==
                                          sign(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ][cville.rank<thr.i & anchor.rank<thr.i]$anchor.beta))
                            st.F <- sum(sign(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ][cville.rank<thr.i & anchor.rank<thr.i]$cville.beta) !=
                                          sign(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ][cville.rank<thr.i & anchor.rank<thr.i]$anchor.beta))
                            
                            bt <- binom.test(st.T, st.T+st.F, .5)
                            
                            
                            tmp1 <- data.table(chr.x=chr.i, 
                                               thr=thr.i , 
                                               perm=k,
                                               anchor.model=model,
                                               win.start=wins$start[pos.ith],
                                               win.end=wins$end[pos.ith],
                                               or=fet$estimate, p=fet$p.value, lci=fet$conf.int[1], uci=fet$conf.int[2],
                                               st.pr=bt$estimate, st.p=bt$p.value, st.lci=bt$conf.int[1], st.uci=bt$conf.int[2])
                            
                            ### return
                            return(tmp1)
                            
                          }
                }
      }
    
    ###save file
    rbind(enrichment.core20, enrichment.cline) -> machado.datasets.enrch
    
    save(machado.datasets.enrch, 
         file = "./machado.datasets.enrch.Rdata"
         # paste("out.enr/",
         #       paste(model, k , "Omnubis.enrich.winlevel", "Rdata",
         #             sep = "."),
         #       sep = "")
    )    
    
    