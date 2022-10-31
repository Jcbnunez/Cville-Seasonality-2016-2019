# ijob -A berglandlab_standard -c1 -p standard --mem=5G
# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

args = commandArgs(trailingOnly=TRUE)
jobId=as.numeric(args[1])

## jobId=1


### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  library(ggplot2)
  library(patchwork)
  library(SeqArray)
  library(metap)
  library(lubridate)

### load inversion
  inv.dt <- fread("~/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr")
  setnames(inv.dt, "stop", "end")
  setkey(inv.dt, start, end)

### common SNP.dt
  load("~/Overwintering_18_19/SNP_filtering/snp_dt_25percMissing.Rdata")
  setkey(snp.dt, id)

### collect results
  perm.sets <- list.files("/scratch/aob2x/dest_glm_morePerms_nested", ".csv", full.name=T)
  header <- perm.sets[grepl("locality", perm.sets)]
  permsets <- perm.sets[!grepl("locality", perm.sets)]

### assign job job
  psi <- perm.sets[jobId]

### samps
  samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
  samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
  samps[nchar(month)==1, month:=paste("0", month, sep="")]
  samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
  samps[nchar(day)==1, day:=paste("0", day, sep="")]
  samps <- samps[set!="dgn"]
  samps[,Date:=date(paste(year, month, day, sep="-"))]

  samps[locality=="UA_od", locality:="UA_Ode"]



  targetLocales <- list("FI_Aka"=list(locality="FI_Aka",   minYear=2000),
                        "DE_Bro"=list(locality="DE_Bro",   minYear=2000),
                        "VA_ch"=list(locality="VA_ch",     minYear=2014),
                        "PA_li"=list(locality="PA_li",     minYear=2000),
                        "DE_Mun"=list(locality="DE_Mun",   minYear=2000),
                        "UA_Ode"=list(locality=c("UA_Ode"), minYear=2010),
                        "TR_Yes"=list(locality="TR_Yes",   minYear=2000))

### load GDS object
  genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)


## get files
  glm.out <- fread(psi)
  glm.header <- fread(header, header=F, nrows=1)
  setnames(glm.out,
            names(glm.out),
            as.character(glm.header[1,]))

  setkey(glm.out, variant.id)

  glm.out <- merge(snp.dt, glm.out, by.x="id", by.y="variant.id")
  glm.out[,p.lrt:= 1 - pchisq(chisq, df)]


### assign ranks: (1) all tested sites; (2) removing rep elements & cM/Mb==0

  lrt.rank.clean <- glm.out[N==0 & cm_mb>0 & !is.na(cm_mb), list(r.clean=rank(-chisq), l.clean=length(chisq),
                                    id),  list(mod, perm)]


  lrt.rank.all <- glm.out[, list(r.all=rank(-chisq), l.all=length(chisq),
                                    id),  list(mod, perm)]
  setkey(lrt.rank.clean, id, mod, perm)
  setkey(lrt.rank.all, id, mod, perm)
  setkey(glm.out, id, mod, perm)

  glm.out <- merge(glm.out, lrt.rank.clean, all.x=T)
  glm.out <- merge(glm.out, lrt.rank.all, all.x=T)

  glm.out[,rnp.clean:=r.clean/(l.clean+1)]
  glm.out[,rnp.all:=r.all/(l.all+1)]

#### FDR
  q.clean <- glm.out[N==0 & cm_mb>0 & !is.na(cm_mb), list(q.clean=p.adjust(p.lrt, "fdr"),
                                    id),  list(mod, perm)]


  q.all <- glm.out[, list(q.all=p.adjust(p.lrt, "fdr"),
                                    id),  list(mod, perm)]
  setkey(q.clean, id, mod, perm)
  setkey(q.all, id, mod, perm)
  setkey(glm.out, id, mod, perm)

  glm.out <- merge(glm.out, q.clean, all.x=T)
  glm.out <- merge(glm.out, q.all, all.x=T)

### save glm.out
#  save(glm.out, file=paste("/scratch/aob2x/summarized_dest_glm/glm.out.", glm.out[1]$locality, "_", glm.out[1]$perm, ".Rdata", sep=""))
  save(glm.out, file=paste("/scratch/aob2x/summarized_dest_glm_nested/glm.out.", glm.out[1]$locality, "_", glm.out[1]$perm, ".Rdata", sep=""))













####################################
### some basic tests and figures ###
####################################

#### temperature test
#  tab.temp.rep <- table(sig=glm.out[mod=="aveTemp"]$rnp.all<0.01,
#              rep=glm.out[mod=="aveTemp"]$N>0)
#  temp.rep <- fisher.test(tab.temp.rep)
#  tab.temp.rec <- table(sig=glm.out[mod=="aveTemp"]$rnp.all<0.01,
#              rec=glm.out[mod=="aveTemp"]$cm_mb==0)
#  temp.rec <- fisher.test(tab.temp.rec)

#### year test
#  tab.year.rep <- table(sig=glm.out[mod=="year_factor"]$rnp.all<0.01, rep=glm.out[mod=="year_factor"]$N>0)
#  year.rep <- fisher.test(tab.year.rep)
#  tab.year.rec <- table(sig=glm.out[mod=="year_factor"]$rnp.all<0.01, rep=glm.out[mod=="year_factor"]$cm_mb==0)
#  year.rec <- fisher.test(tab.year.rec)

#### all inversion, all tests
#  fet.invs <- foreach(inv.i=c(inv.dt$invName), .combine="rbind")%dopar%{
#    # inv.i="2Lt"
#    foreach(mod.i=c("aveTemp", "year_factor"), .combine="rbind")%do%{
#      print(paste(inv.i, mod.i, sep=" / "))
#      # mod.i="aveTemp"
#      tab <- table(sig=glm.out[mod==mod.i]$rnp.clean<0.01,
#                  Inv=grepl(inv.i, glm.out[mod==mod.i]$invName))

#      temp.inv <- fisher.test(tab)
#      data.table(model=mod.i, factor=inv.i,
#                or=temp.inv$estimate, or.lci=temp.inv$conf.int[1], or.uci=temp.inv$conf.int[2],
#                FF=tab[1,1], FT=tab[1,2], TF=tab[2,1], TT=tab[2,2])
#    }
#  }

#### merge
#  fet.dt <- rbind(
#            fet.invs,
#            data.table(model=c(rep("aveTemp", 2), rep("year_factor", 2)),
#                       factor=rep(c("Repetitive", "LowRec"), 2),
#                        or=c(temp.rep$estimate, temp.rec$estimate,
#                             year.rep$estimate, year.rec$estimate),
#                        or.lci=c(temp.rep$conf.int[1], temp.rec$conf.int[1],
#                                year.rep$conf.int[1], year.rec$conf.int[1]),
#                        or.uci=c(temp.rep$conf.int[2], temp.rec$conf.int[2],
#                                 year.rep$conf.int[2], year.rec$conf.int[2]),
#                        FF=c(tab.temp.rep[1,1], tab.temp.rec[1,1],
#                             tab.year.rep[1,1], tab.year.rec[1,1]),
#                        FT=c(tab.temp.rep[1,2], tab.temp.rec[1,2],
#                             tab.year.rep[1,2], tab.year.rec[1,2]),
#                        TF=c(tab.temp.rep[2,1], tab.temp.rec[2,1],
#                             tab.year.rep[2,1], tab.year.rec[2,1]),
#                        TT =c(tab.temp.rep[2,2], tab.temp.rec[2,2],
#                              tab.year.rep[2,2], tab.year.rec[2,2]))
#            )

#  fet.dt[,psi:=psi]


#### Q. fet
#fet.invs <- foreach(inv.i=c(inv.dt$invName), .combine="rbind", .errorhandling="remove")%dopar%{
#  # inv.i="2Lt"
#  foreach(mod.i=c("aveTemp", "year_factor"), .combine="rbind", .errorhandling="remove")%do%{
#    print(paste(inv.i, mod.i, sep=" / "))
#    # mod.i="aveTemp"
#    tab <- table(sig=glm.out[mod==mod.i]$q.clean<=0.1,
#                Inv=grepl(inv.i, glm.out[mod==mod.i]$invName))

#    temp.inv <- fisher.test(tab)
#    data.table(model=mod.i, factor=inv.i,
#              or=temp.inv$estimate, or.lci=temp.inv$conf.int[1], or.uci=temp.inv$conf.int[2], p=temp.inv$p.value,
#              FF=tab[1,1], FT=tab[1,2], TF=tab[2,1], TT=tab[2,2])
#  }
#}
#fet.invs[,psi:=psi]

# fet.invs <- foreach(inv.i=c(inv.dt$invName), .combine="rbind", .errorhandling="remove")%do%{
#   # inv.i="2Lt"
#   foreach(rnp=expand.grid(sapply(c(1:9), function(x) x*10^c(-4:-1)))[,1], .combine="rbind", .errorhandling="remove")%dopar%{
#     print(paste(inv.i, rnp, sep=" / "))
#     mod.i="aveTemp"
#     tab <- table(sig=glm.out[mod==mod.i]$rnp.clean<=rnp,
#                 Inv=grepl(inv.i, glm.out[mod==mod.i]$invName))
#
#     temp.inv <- fisher.test(tab)
#     data.table(model=mod.i, factor=inv.i, i=rnp, psi=psi,
#               or=temp.inv$estimate, or.lci=temp.inv$conf.int[1], or.uci=temp.inv$conf.int[2],
#               p=temp.inv$p.value,
#               FF=tab[1,1], FT=tab[1,2], TF=tab[2,1], TT=tab[2,2])
#   }
# }
# fet.invs




# #############
# # windows ###
# #############

#  define windows
# win.bp <- 1e5
# step.bp <- 5e4
# setkey(snp.dt, "chr")
# wins <- foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind", .errorhandling="remove")%do%{
#     tmp <- snp.dt[J(chr.i)]
#     data.table(chr=chr.i,
#                 start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
#                 end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
# }
# wins[,i:=1:dim(wins)[1]]
# dim(wins)

#  run windows
# setkey(glm.out, chr, pos)

# win.out <- foreach(win.i=c(1:dim(wins)[1]), .errorhandling="remove")%dopar%{
#   # win.i <- 1
#   message(paste(win.i, dim(wins)[1], sep=" / "))
#   win.tmp <- glm.out[J(data.table(chr=wins[win.i]$chr, pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]
#   #win.tmp <- glm.out[J(data.table(index=wins[win.i]$start:wins[win.i]$end, key="index")), nomatch=0]
#   win.tmp[,Z:=qnorm(p.lrt, 0, 1)]
#   win.tmp[,rnpZ:=qnorm(rnp.clean, 0, 1)]

#   seqSetFilter(genofile, variant.id=unique(win.tmp$id),
#             sample.id=samps[
#                             locality==targetLocales[which(psi==names(targetLocales))][[1]]$locality
#                             ][
#                             year>=targetLocales[which(psi==names(targetLocales))][[1]]$minYear]$sampleId)

#   af <- seqGetData(genofile, "annotation/format/FREQ")
#   f.hat <- data.table(fhat=colMeans(af[[2]], na.rm=T), id=seqGetData(genofile, "variant.id"))

#   win.tmp <- merge(win.tmp, f.hat, by="id")
#   win.tmp[,het:=fhat*(1-fhat)]

#   win.tmp[!is.na(rnp.clean),
#             list(rnp.pr=mean(rnp.clean<=0.05),
#                 pr=c(0.05),
#                 rbinom.p=binom.test(sum(rnp.clean<=0.05), length(rnp.clean), .05)$p.value,
#                 wZa=sum(het*Z)/(sqrt(sum(het^2))),
#                 wZa.p=pnorm(sum(het*Z)/(sqrt(sum(het^2))), lower.tail=T),
#                 rnp.wZa=sum(het*rnpZ)/(sqrt(sum(het^2))),
#                 rnp.wZa.p=pnorm(sum(het*rnpZ)/(sqrt(sum(het^2))), lower.tail=T),
#                 psi=psi,
#                 mean.chisq=c(mean(chisq)),
#                 median.chisq=c(median(chisq)),
#                 min.p.lrt=min(p.lrt),
#                 min.rank=min(r.clean),
#                 min.rnp=min(rnp.clean),
#                 gr=c("obs"),
#                 nSNPs=.N, i=win.i,
#                 start.bp=min(win.tmp$pos), stop.bp=max(win.tmp$pos), win.i=win.i),
#             list(mod, chr, invName, delta, locality)]
# }
# win.out <- rbindlist(win.out)
# win.out <- merge(win.out, wins, by="i")
# win.out[mod=="aveTemp"][order(Z=wZa)]
# win.out[mod=="aveTemp"][order(rbinom.p)]

### save

  #save(fet.dt, file=paste("/scratch/aob2x/summarized_dest_glm/fet.dt.", psi, ".Rdata", sep=""))
  #save(win.out, win.out.exp, file=paste("/scratch/aob2x/summarized_dest_glm/window.", psi, ".Rdata", sep=""))
  #save(fet.invs, file=paste("/scratch/aob2x/summarized_dest_glm/fetInvs.qb.singleLocale.dt.", psi, ".Rdata", sep=""))
  #save(win.out, file=paste("/scratch/aob2x/summarized_dest_glm/WZA_window.dt.", psi, ".Rdata", sep=""))







  #### run windows
  #  setkey(glm.out, chr, pos)
  #  win.out <- foreach(win.i=c(1:dim(wins)[1]), .errorhandling="remove")%dopar%{
  #    # win.i <- 126
  #    message(paste(win.i, dim(wins)[1], sep=" / "))
  #    win.tmp <- glm.out[J(data.table(chr=wins[win.i]$chr, pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]
  #
  #
  #    win.tmp[!is.na(rnp.clean),
  #              list(rnp.pr=c(mean(rnp.clean<=0.05), mean(rnp.clean<=0.01)),
  #                  fm.p= 1 - c(p.tfisher(q=stat.tfisher(p=p.lrt, tau1=0.05, tau2=0.25), n=length(p.lrt), tau1=0.05, tau2=.25),
  #                              p.tfisher(q=stat.tfisher(p=p.lrt, tau1=0.01, tau2=0.25), n=length(p.lrt), tau1=0.01, tau2=.25)),
  #                  fm.stat= c(stat.tfisher(p=p.lrt, tau1=0.05, tau2=.25),
  #                         stat.tfisher(p=p.lrt, tau1=0.01, tau2=.25)),
  #                  rtp=c(ranktruncated(p.lrt, K=25)$RTP$p.Value) ,
  #                  pr=c(0.05, .01),
  #                  mean.chisq=c(mean(chisq)),
  #                  median.chisq=c(median(chisq)),
  #                  gr=c("obs"),
  #                  nSNPs=.N, i=win.i),
  #              list(mod, chr, invName, delta, locality)]
  #
  #  }
  #  win.out <- rbindlist(win.out)
  #  win.out <- merge(win.out, wins, by="i")
  #  win.out[mod=="aveTemp"][order(rtp)]
  #
  #
  #  win.out.exp <- data.table(median.chisq=qchisq(.5, c(1,4)),
  #                                uci.chisq=qchisq(.025, c(1,4)),
  #                                lci.chisq=qchisq(.975, c(1,4)),
  #                                mod=c("aveTemp", "year_factor"))
