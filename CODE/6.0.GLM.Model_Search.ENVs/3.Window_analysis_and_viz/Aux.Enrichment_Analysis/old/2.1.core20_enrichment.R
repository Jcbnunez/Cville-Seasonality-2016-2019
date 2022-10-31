# ijob -A berglandlab_standard -c1 -p standard --mem=8G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

  args = commandArgs(trailingOnly=TRUE)
  #job=as.numeric(args[1])-1
  z=as.numeric(args[1])
  
### libraries
  library(data.table)
  library(tidyverse)
  library(magrittr)
  library(foreach)
  library(doMC)
  registerDoMC(4)
  library(tidyr)

### load core20 GLM object
  message("load CORE 20 GLM")
  load(file="/project/berglandlab/alan/core20glm.Rdata")

### load thermal GLM object
  sets <- data.table(mod=c(1:11),
                     start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                     end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))
  
#job <- 0
##  glm.fn <- paste("/project/berglandlab/thermal_glm_dest/dest_glm_final_nested_qb/processedGLM##/glm.out.VA_ch_", job, ".Rdata", sep="")
##  load(glm.fn)

  ### load VA GLM object
  #model <- "temp.max_2"
  #variable <- "temp.max"
  #mod.num <- 2 
  #perm_ith <- 0
  
  ### guide_file
  mods_guide = 
  data.frame(
    models= c("temp.max_2", "temp.ave_4", "precip.ave_11", "precip.var_3"),
    variables= c("temp.max","temp.ave", "precip.ave", "precip.var"),
    mod.num= c(2,4,11,3)
      )
  
  perm_ith <- 0
  
    model <- mods_guide$models[z]
    variable <- mods_guide$variables[z]
    mod.num <-  mods_guide$mod.num[z]
    perm_ith <-  perm_ith
  
  message("Ingestigating the following model=")
  message(paste(model, variable,mod.num,  perm_ith, sep = " / "))  
    
  glm.fn <- paste("/project/berglandlab/alan/environmental_ombibus/",
                  model, 
                  "/",
                  paste(model, "glmRNP", sep = "."), 
                  ".Rdata", 
                  sep="")
  
  message("load Cville GLM for selected model")
  load(glm.fn)
  ## --> 
  glm.out.0 <- glm.out %>%
                filter(perm == perm_ith)
  
### merge
  message("Merge datasets")
  
  m1 <- merge(glm.out.0[variable==variable], 
              core20.glm[mod=="season+locality_factor"][set=="machado"], 
              by.x="variant.id", by.y="variant.id")
  
  m2 <- merge(glm.out.0[variable==variable], 
              core20.glm[mod=="season+locality_factor"][set=="no_va"],
              by.x="variant.id", by.y="variant.id")

### clean up
  rm(glm.out, core20.glm)

### merge with haplotags
#  load("/project/berglandlab/jcbnunez/Shared_w_Alan/haplo_tags_SNPids.Rdata")
#  haplo_tags_SNPids_and_inv <- as.data.table(haplo_tags_SNPids_and_inv)
#  haplo_tags_SNPids_and_inv[,pos:=as.numeric(pos)]
#  setkey(m1, chr, pos)
#  setkey(m2, chr, pos)
#
#  setkey(haplo_tags_SNPids_and_inv, chr, pos)
#
#  m1 <- merge(m1, haplo_tags_SNPids_and_inv, all.x=T)
#  m2 <- merge(m2, haplo_tags_SNPids_and_inv, all.x=T)

### rank normalization
  m1 <- m1[chr!="X"]
  m2 <- m2[chr!="X"]

  # x --> thermal is p_lrt
  # y --> seasonal is p.lrt
  
  m1[!is.na(p_lrt) & !is.na(p.lrt) ,thermal.rank := rank(p_lrt)/length(p_lrt)]
  m1[!is.na(p_lrt) & !is.na(p.lrt) ,season.rank := rank(p.lrt)/length(p.lrt)]


  m2[!is.na(p_lrt) & !is.na(p.lrt) ,thermal.rank := rank(p_lrt)/length(p_lrt)]
  m2[!is.na(p_lrt) & !is.na(p.lrt) ,season.rank := rank(p.lrt)/length(p.lrt)]


### get betas
 # m1.beta <- m1[,list(thermal.beta=tstrsplit(b.x, ";")%>%
 #                       last%>%
 #                       as.numeric, 
 #                     seas.beta=tstrsplit(b.y, ";")%>%last%>%as.numeric), 
 #               list(id)]
  message("Estimating the betas")
  
  m1.beta <- m1[,list(thermal.beta=b_temp%>%
                        last%>%
                        as.numeric, 
                      seas.beta=tstrsplit(b, ";")%>%
                        last%>%
                        as.numeric), 
                list(variant.id)]
  
  
  m2.beta <- m2[,list(thermal.beta=b_temp%>%
                        last%>%
                        as.numeric,
                      seas.beta=tstrsplit(b, ";")%>%
                        last%>%
                        as.numeric),
                list(variant.id)]

  
  m1 <- merge(m1, m1.beta, by="variant.id")
  m2 <- merge(m2, m2.beta, by="variant.id")


### enrichment test
  m1[,inv:=invName!="none"]
  m2[,inv:=invName!="none"]

### save
  #save(m1, m2, file="/scratch/aob2x/core20_dest_aveTemp.Rdata")
  #thrs <- expand.grid(sapply(1:9, function(x) x*10^(-5:-1)))[,1]
  thrs <- c(
            #0.01, 
            0.05)
  
  message("Running Enrichment")
  
  o <- foreach(chr.i=unique(m1$chr), 
               .combine="rbind", 
               .errorhandling="pass")%do%{
                 
      foreach(inv.i=c(T,F), 
              .combine="rbind", 
              .errorhandling="pass")%do%{
                
      foreach(thr.i=thrs, 
              .combine="rbind", 
              .errorhandling="pass")%do%{
        
        # chr.i <- "3R"; inv.i <- T; thr.i<-0.05
          message(paste(chr.i, inv.i, thr.i, sep=" / "))

          ### Machado set, enrichment
            tab <- table(m1[chr==chr.i][inv==inv.i]$thermal.rank  < thr.i,
                         m1[chr==chr.i][inv==inv.i]$season.rank <  thr.i)

            fet <- fisher.test(tab)

         ### Machado set, sign test
           st.T <- sum(sign(m1[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$thermal.beta) ==
                        sign(m1[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$seas.beta))
           st.F <- sum(sign(m1[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$thermal.beta) !=
                        sign(m1[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$seas.beta))

           bt <- binom.test(st.T, st.T+st.F, .5)


           tmp1 <- data.table(chr=chr.i, inv=inv.i, thr=thr.i, perm=perm_ith,
                      mod="core20_Machado", set=model,
                      or=fet$estimate, p=fet$p.value, lci=fet$conf.int[1], uci=fet$conf.int[2],
                      st.pr=bt$estimate, st.p=bt$p.value, st.lci=bt$conf.int[1], st.uci=bt$conf.int[2])


         ### DEST set, enrichment
           tab <- table(m2[chr==chr.i][inv==inv.i]$thermal.rank  < thr.i,
                        m2[chr==chr.i][inv==inv.i]$season.rank <  thr.i)

           fet <- fisher.test(tab)

         ### DEST set, sign test
           st.T <- sum(sign( m2[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$thermal.beta) ==
                        sign(m2[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$seas.beta))
           st.F <- sum(sign( m2[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$thermal.beta) !=
                        sign(m2[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$seas.beta))

           bt <- binom.test(st.T, st.T+st.F, .5)


           tmp2 <- data.table(chr=chr.i, inv=inv.i, thr=thr.i, perm=perm_ith,
                      mod="core20_DEST", set=model,
                      or=fet$estimate, p=fet$p.value, lci=fet$conf.int[1], uci=fet$conf.int[2],
                      st.pr=bt$estimate, st.p=bt$p.value, st.lci=bt$conf.int[1], st.uci=bt$conf.int[2])

        ### return
          rbind(tmp1, tmp2, fill=T)

      }
    }
  }


  message("Save File")
  
save(o, 
     file = 
       
       paste(model, "core20.enrrich", "Rdata",
             sep = ".")
       
       )    
  
#  save(o.dest_core20, file=paste("/project/berglandlab/alan/core20_enrichment/dest_core20.enrichment_signTest.", job, ".Rdata"))

  
#### WINDOW - haplotag SNPs
#  m1[season.rank<.05,list(pr=mean(sign(thermal.beta)==sign(seas.beta.y), na.rm=T), .N), list(win)]
#  m2[season.rank<.05,list(pr=mean(sign(thermal.beta)==sign(seas.beta), na.rm=T), .N), list(win)]
#
#  m2[,list(TT=sum(season.rank<.05 & thermal.rank<.05),
#          TF=sum(season.rank<.05 & thermal.rank>.05),
#           FT=sum(season.rank>.05 & thermal.rank<.05),
#          FF=sum(season.rank>.05 & thermal.rank>.05)), list(win)]
##
##  o.dest_core20 %>%
##    ggplot(aes(
##      x=chr,
##      y=or,
##      ymin=lci,
##      ymax=uci,
##      #color = mod,
##      shape = inv
##    )) +
##    geom_hline(yintercept =1) +
##    geom_point(size = 4, position=position_dodge(width=0.5)) +
##    scale_y_continuous(trans='log2') +
##    geom_errorbar(
##      width = 0.5,
##      linetype = "dotted",
##      position=position_dodge(width=0.5)) +
##    facet_grid(~mod) -> core_20_temp2
##  
##  ggsave(core_20_temp2, file = "core_20_Tave4.dest.pdf", h = 3)
##  
##
#
#  m1[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,thermal.rank := rank(p.lrt.x)/length(p.lrt.x)]
#  m1[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,season.rank := rank(p.lrt.y)/length(p.lrt.y)]
#
#
#  m2[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,thermal.rank := rank(p.lrt.x)/length(p.lrt.x)]
#  m2[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,season.rank := rank(p.lrt.y)/length(p.lrt.y)]
