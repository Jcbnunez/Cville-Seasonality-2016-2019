# ijob -A berglandlab_standard -c5 -p standard --mem=20G
# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

args = commandArgs(trailingOnly=TRUE)
jobId=as.numeric(args[1])

## jobId=67

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)

### collect results
  perm.sets <- list.files("/scratch/aob2x/dest_glm", full.name=F)

### assign job job
  psi <- perm.sets[jobId]

### thrs
  thrs <- expand.grid(sapply(c(1:9), function(x) x*10^-c(3:1)))[,1]
  #thrs <- .05

### load file
  load(file=paste("/scratch/aob2x/summarized_dest_glm/glm.out.", psi, ".Rdata", sep=""))
  glm.out <- glm.out[N==0 & cm_mb>0 & !is.na(cm_mb)]

    #q.clean <- glm.out[N==0 & cm_mb>0 & !is.na(cm_mb), list(q.clean=p.adjust(p.lrt, "fdr"),
  #                                  id),  list(mod)]
#
#
  #q.all <- glm.out[, list(q.all=p.adjust(p.lrt, "fdr"),
  #                                  id),  list(mod)]
  #setkey(q.clean, id, mod)
  #setkey(q.all, id, mod)
  #setkey(glm.out, id, mod)
#
  #glm.out <- merge(glm.out, q.clean, all.x=T)
  #glm.out <- merge(glm.out, q.all, all.x=T)

### load eLife seasonal SNps
  elife <- fread("/project/berglandlab/alan/drosRTEC/mel_all_paired20_2sample_caF_popyear.f_s.dm6.bed", header=F)
  setnames(elife, names(elife), c("chr", "pos", "V3", "b", "p", "V6"))
  elife[,chr:=gsub("chr", "", chr)]

  setkey(elife, chr, pos)
  setkey(glm.out, chr, pos)

### merge
  m <- merge(glm.out[mod=="aveTemp"], elife)
  m[,temp.rank:=rank(p.lrt)/(length(p.lrt)+1)]
  m[,elife.rank:=rank(p.y)/(length(p.y)+1)]

  setkey(m, chr)
  elife.ft <- foreach(i=thrs, .combine="rbind")%dopar%{

    tab <- table(m$temp.rank<=i, m$elife.rank<=i)
    ft <- fisher.test(tab)
    dt <- data.table(i=i, or=ft$estimate, p=ft$p.value, lci=ft$conf.int[1], uci=ft$conf.int[2], chr="genome",
                    FF=tab[1,1], FT=tab[1,2], TF=tab[2,1], TT=tab[2,2])

    dtc <-foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind")%do%{
      message(paste(i, chr.i, sep=" / "))
      m.tmp <- m[J(chr.i)]
      tab <- table(m.tmp$temp.rank<=i, m.tmp$elife.rank<=i)
      ft <- fisher.test(tab)
      data.table(i=i, or=ft$estimate, p=ft$p.value, lci=ft$conf.int[1], uci=ft$conf.int[2], chr=chr.i,
                      FF=tab[1,1], FT=tab[1,2], TF=tab[2,1], TT=tab[2,2])

    }
    rbind(dt, dtc)
  }
  elife.ft[order(i)]
  elife.ft[,paper:="eLife"]

### 2014 data
  #load("/project/berglandlab/alan/drosRTEC/seas2014.Rdata")
  #seas2014 <- as.data.table(seas2014)
  #write.table(seas2014, quote=F, row.names=F, col.names=F, sep="\t", file="/project/berglandlab/alan/drosRTEC/seas2014.delim")

  # Overwintering_18_19/liftOver_PlosG_eLife_seasonalSNPs/liftOver_plosG.sh

  pg <- fread("/project/berglandlab/alan/drosRTEC/seas2014.dm6.bed")
  setnames(pg, names(pg), c("chr", "pos", "V3", "sq"))
  pg[,chr:=gsub("chr", "", chr)]

  setkey(pg, chr, pos)
  setkey(glm.out, chr, pos)

### merge
  m2 <- merge(glm.out[mod=="aveTemp"], pg)
  m2[,temp.rank:=rank(p.lrt)/(length(p.lrt)+1)]
  m2[,pg.rank:=rank(sq)/(length(sq)+1)]

  setkey(m2, chr)
  pg.ft <- foreach(i=thrs, .combine="rbind")%dopar%{

    tab <- table(m2$temp.rank<=i, m2$pg.rank<=i)
    ft <- fisher.test(tab)
    dt <- data.table(i=i, or=ft$estimate, p=ft$p.value, lci=ft$conf.int[1], uci=ft$conf.int[2], chr="genome",
                    FF=tab[1,1], FT=tab[1,2], TF=tab[2,1], TT=tab[2,2])

    dtc <-foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind")%do%{
      message(paste(i, chr.i, sep=" / "))
      m.tmp <- m2[J(chr.i)]
      tab <- table(m.tmp$temp.rank<=i, m.tmp$pg.rank<=i)
      ft <- fisher.test(tab)
      data.table(i=i, or=ft$estimate, p=ft$p.value, lci=ft$conf.int[1], uci=ft$conf.int[2], chr=chr.i,
                      FF=tab[1,1], FT=tab[1,2], TF=tab[2,1], TT=tab[2,2])

    }
    rbind(dt, dtc)
  }
  pg.ft[order(i)]
  pg.ft[,paper:="PlosG"]

### combine
  o <- rbind(elife.ft, pg.ft)
  o[,psi:=psi]

### save
  #system("mkdir /scratch/aob2x/temp_seasonal_SNP_overlap")
  save(o, file=paste("/scratch/aob2x/temp_seasonal_SNP_overlap/", psi, ".Rdata", sep=""))
