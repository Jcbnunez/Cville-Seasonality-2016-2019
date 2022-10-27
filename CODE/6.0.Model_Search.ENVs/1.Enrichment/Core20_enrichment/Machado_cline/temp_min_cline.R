# ijob -A berglandlab_standard -c5 -p standard --mem=40G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)

### load VA GLM object
  job <- 0
  glm.fn <- paste("/project/berglandlab/alan/environmental_ombibus/temp.ave_4/temp.ave_4.glmRNP.Rdata", sep="")
  load(glm.fn)

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

  core20 <- merge(core20, liftover)

  setnames(core20, c("dm6_chr", "dm6_pos"), c("chr", "pos"))

### merge with VA_glm
  setkey(core20, chr, pos)
  setkey(glm.out, chr, pos)
  m <- merge(glm.out, core20)

  rm(glm.out, core20, liftover)

### load premade
  # load("/scratch/aob2x/drosRTEC_DEST_merge.Rdata")

### quick test
  setkey(m, perm)
  o <- foreach(perm.i=0:100)%do%{
      #perm.i <- 0
      message(perm.i)
      tmp <- m[J(perm.i)][!is.na(rnp) & !is.na(clinal.p)]

      tmp[,rnp.clean:=rank(p_lrt)/(length(p_lrt+1))]
      tmp[,rnp.machado:=rank(clinal.p)/(length(clinal.p+1))]
      tmp[,inv:=invName!="none"]


      #thrs <- expand.grid(sapply(1:9, function(x) x*10^(-5:-1)))[,1]
      thrs <- c(0.01, 0.05)

      o <- foreach(chr.i=unique(tmp$chr), .combine="rbind", .errorhandling="pass")%do%{
        foreach(inv.i=c(T,F), .combine="rbind", .errorhandling="pass")%do%{
          foreach(thr.i=thrs, .combine="rbind", .errorhandling="pass")%dopar%{
              # chr.i <- "2R"; inv.i=T; thr.i=0.01
              #message(paste(mod.i, chr.i, inv.i, thr.i, sep=" / "))
              tab <- table(tmp[chr==chr.i][inv==inv.i]$rnp.clean  <thr.i,
                           tmp[chr==chr.i][inv==inv.i]$rnp.machado<thr.i)
              fet <- fisher.test(tab)


              st.T <- sum(sign( tmp[chr==chr.i][inv==inv.i][rnp.clean<thr.i & rnp.machado<thr.i]$b_temp) ==
                           sign(tmp[chr==chr.i][inv==inv.i][rnp.clean<thr.i & rnp.machado<thr.i]$clinal.coef))
              st.F <- sum(sign( tmp[chr==chr.i][inv==inv.i][rnp.clean<thr.i & rnp.machado<thr.i]$b_temp) !=
                           sign(tmp[chr==chr.i][inv==inv.i][rnp.clean<thr.i & rnp.machado<thr.i]$clinal.coef))

              bt <- binom.test(st.T, st.T+st.F, .5)


              data.table(chr=chr.i, inv=inv.i, thr=thr.i,
                         perm.i=perm.i, or=fet$estimate, p=fet$p.value, lci=fet$conf.int[1], uci=fet$conf.int[2],
                          st.pr=bt$estimate, st.p=bt$p.value, st.lci=bt$conf.int[1], st.uci=bt$conf.int[2], nTT=st.T+st.F)
          }
        }
      }
      return(o)
    }

o <- rbindlist(o)

o.cline <- o

#### save
  save(o.cline, file="~/drosRtec_enrichment_cline.Rdata")

### scp
  scp aob2x@rivanna.hpc.virginia.edu:~/drosRtec_enrichment.Rdata ~/.

### libraries
  library(ggplot2)
  library(data.table)

  load("~/drosRtec_enrichment.Rdata")

  ggplot() +
  geom_line(data= o, aes(x=log10(thr), y=or, group=interaction(mod, set), color=mod, linetype=set)) +
  geom_point(data=o[p<.0005], aes(x=log10(thr), y=or, group=interaction(mod, set), color=mod, linetype=set)) +
  #geom_point(data=o[p<.005], aes(x=interaction(mod, set), y=log2(or), group=thr), color="black", size=.5) +
  facet_grid(chr~inv) +
  theme(axis.text.x=element_text(angle=90))
