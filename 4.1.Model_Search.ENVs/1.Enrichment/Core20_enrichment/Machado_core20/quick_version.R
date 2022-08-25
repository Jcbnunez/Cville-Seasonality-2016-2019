# ijob -A berglandlab_standard -c5 -p standard --mem=20G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)

### load VA GLM object
  job <- 0
  glm.fn <- paste("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_", job, ".Rdata", sep="")
  load(glm.fn)

### load core20 data
  core20.orig <- fread("/project/berglandlab/alan/drosRTEC/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear.f_s.glm")
  core20.swap <- fread("/project/berglandlab/alan/drosRTEC/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm")

  core20.orig[,set:="orig"]
  core20.swap[,set:="swap"]
  core20 <- rbind(core20.orig, core20.swap)

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

### load premade
  load("/scratch/aob2x/drosRTEC_DEST_merge.Rdata")

### quick test
  setkey(m, set, mod)
  o <- foreach(set.i=c("aveTemp+year_factor", "year_factor"), .combine="rbind", .errorhandling="remove")%do%{
    foreach(mod.i=c("orig", "swap"), .combine="rbind", .errorhandling="remove")%do%{
      #mod.i <- "swap"; set.i<-"year_factor"
      tmp <- m[J(data.table(set=set.i, mod=mod.i, key="mod,set"))][!is.na(rnp.clean) & !is.na(seas.p)]

      tmp[,rnp.clean:=rank(p.lrt)/(length(p.lrt+1))]
      tmp[,rnp.machado:=rank(seas.p)/(length(seas.p+1))]
      tmp[,inv:=invName!="none"]
      thrs <- expand.grid(sapply(1:9, function(x) x*10^(-5:-1)))[,1]
      #thrs <- c(0.01, 0.05)

      o <- foreach(chr.i=unique(tmp$chr), .combine="rbind", .errorhandling="remove")%do%{
        foreach(inv.i=c(T,F), .combine="rbind", .errorhandling="remove")%do%{
          foreach(thr.i=thrs, .combine="rbind", .errorhandling="remove")%dopar%{
              message(paste(set.i, mod.i, chr.i, inv.i, thr.i, sep=" / "))
              tab <- table(tmp[chr==chr.i][inv==inv.i]$rnp.clean  <thr.i,
                           tmp[chr==chr.i][inv==inv.i]$rnp.machado<thr.i)
              fet <- fisher.test(tab)

              data.table(chr=chr.i, inv=inv.i, thr=thr.i,
                         mod=mod.i, set=set.i, or=fet$estimate, p=fet$p.value, lci=fet$conf.int[1], uci=fet$conf.int[2])
          }
        }
      }
      return(o)
    }
}


#### save
  save(o, file="~/drosRtec_enrichment.Rdata")

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
