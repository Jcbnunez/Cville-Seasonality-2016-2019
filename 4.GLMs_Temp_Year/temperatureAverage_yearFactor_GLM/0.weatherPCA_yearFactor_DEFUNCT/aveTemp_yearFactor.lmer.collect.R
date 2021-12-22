# ijob -A berglandlab_standard -c5 -p standard --mem=20G
# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

args = commandArgs(trailingOnly=TRUE)
jobId=as.numeric(args[1])

## jobId=2


### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  library(ggplot2)
  library(patchwork)


### load inversion
  inv.dt <- fread("~/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr")
  setnames(inv.dt, "stop", "end")
  setkey(inv.dt, start, end)

### common SNP.dt
  load("~/Overwintering_18_19/SNP_filtering/snp_dt_25percMissing.Rdata")
  setkey(snp.dt, id)

### collect results
  perm.sets <- list.files("/scratch/aob2x/dest_glmer", full.name=F)

### assign job job
  psi <- perm.sets[jobId]

### get files
  fl <- list.files(paste("/scratch/aob2x/dest_glmer/", psi, sep=""), full.name=T)
  glm.out <- foreach(fl.i=fl)%dopar%{
    print(fl.i)
    print(fl.i)
    #fl.i <- fl[1]
    load(fl.i)
    return(tmp)
  }
  glm.out <- rbindlist(glm.out)
  setkey(glm.out, variant.id)

  glm.out <- merge(snp.dt, glm.out, by.x="id", by.y="variant.id")
  #glm.out[,p.lrt:= 1 - pchisq(chisq, df)]


### assign ranks: (1) all tested sites; (2) removing rep elements & cM/Mb==0

  lrt.rank.clean <- glm.out[N==0 & cm_mb>0 & !is.na(cm_mb) & nSamps>50, list(r.clean=rank(-chisq), l.clean=length(chisq),
                                    id),  list(mod)]


  lrt.rank.all <- glm.out[, list(r.all=rank(-chisq), l.all=length(chisq),
                                    id),  list(mod)]
  setkey(lrt.rank.clean, id, mod)
  setkey(lrt.rank.all, id, mod)
  setkey(glm.out, id, mod)

  glm.out <- merge(glm.out, lrt.rank.clean, all.x=T)
  glm.out <- merge(glm.out, lrt.rank.all, all.x=T)

  glm.out[,rnp.clean:=r.clean/(l.clean+1)]
  glm.out[,rnp.all:=r.all/(l.all+1)]

#### load precomputed files


  q.clean <- glm.out[N==0 & cm_mb>0 & !is.na(cm_mb) & nSamps>50, list(q.clean=p.adjust(p.lrt, "fdr"),
                                    id),  list(mod)]


  q.all <- glm.out[, list(q.all=p.adjust(p.lrt, "fdr"),
                                    id),  list(mod)]
  setkey(q.clean, id, mod)
  setkey(q.all, id, mod)
  setkey(glm.out, id, mod)

  glm.out <- merge(glm.out, q.clean, all.x=T)
  glm.out <- merge(glm.out, q.all, all.x=T)


####################################
### some basic tests and figures ###
####################################
  fet.invs <- foreach(inv.i=c(inv.dt$invName), .combine="rbind", .errorhandling="remove")%do%{
    # inv.i="2Lt"
    foreach(rnp=expand.grid(sapply(c(1:9), function(x) x*10^c(-4:-1)))[,1], .combine="rbind", .errorhandling="remove")%dopar%{
      print(paste(inv.i, rnp, sep=" / "))
      mod.i="aveTemp"
      tab <- table(sig=glm.out[mod==mod.i]$rnp.clean<=rnp,
                  Inv=grepl(inv.i, glm.out[mod==mod.i]$invName))

      temp.inv <- fisher.test(tab)
      data.table(model=mod.i, factor=inv.i, i=rnp,
                or=temp.inv$estimate, or.lci=temp.inv$conf.int[1], or.uci=temp.inv$conf.int[2],
                p=temp.inv$p.value,
                FF=tab[1,1], FT=tab[1,2], TF=tab[2,1], TT=tab[2,2])
    }
  }
  fet.invs
  #fet.invs[factor=="2Lt"][order(i)]


###############
### windows ###
###############
  #### define windows
    win.bp <- 5e4
    step.bp <- 1e4
    setkey(snp.dt, "chr")
    wins <- foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind", .errorhandling="remove")%do%{
        tmp <- snp.dt[J(chr.i)]
        data.table(chr=chr.i,
                    start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                    end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
    }
    wins[,i:=1:dim(wins)[1]]
    dim(wins)

### run windows
  setkey(glm.out, chr, pos)
  win.out <- foreach(win.i=c(1:dim(wins)[1]), .errorhandling="remove")%dopar%{
    # win.i <- 126
    message(paste(win.i, dim(wins)[1], sep=" / "))
    win.tmp <- glm.out[J(data.table(chr=wins[win.i]$chr, pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]


    win.tmp[!is.na(rnp.clean),
              list(rnp.pr=c(mean(rnp.clean<=0.05), mean(rnp.clean<=0.01)),
                   binom.p=c(binom.test(sum(rnp.clean<=0.05), .N, 0.05)$p.value,
                              binom.test(sum(rnp.clean<=0.01), .N, 0.01)$p.value),
                    pr=c(0.05, .01),
                  mean.chisq=c(mean(chisq)),
                  median.chisq=c(median(chisq)),
                  gr=c("obs"),
                  nSNPs=.N, i=win.i),
              list(mod, chr, invName, delta, locality)]

  }
  win.out <- rbindlist(win.out)
  win.out <- merge(win.out, wins, by="i")
  win.out[mod=="aveTemp"][order(binom.p)]

### save
  #save(glm.out, file=paste("/scratch/aob2x/summarized_dest_glm/glm.out.", psi, ".Rdata", sep=""))
  #save(fet.dt, file=paste("/scratch/aob2x/summarized_dest_glm/fet.dt.", psi, ".Rdata", sep=""))
  #save(win.out, win.out.exp, file=paste("/scratch/aob2x/summarized_dest_glm/window.", psi, ".Rdata", sep=""))
  save(fet.invs, glm.out, win.out, file=paste("/scratch/aob2x/summarized_dest_glm/lmer.all.", psi, ".Rdata", sep=""))
  save(fet.invs, file=paste("/scratch/aob2x/summarized_dest_glm/lmer.fetInv.", psi, ".Rdata", sep=""))
