# ijob -A berglandlab_standard -c5 -p standard --mem=20G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


args = commandArgs(trailingOnly=TRUE)
jobId=as.numeric(args[1])

## jobId=607
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

### load glm.out
  fl <- list.files("/scratch/aob2x/dest_glm_final_nested_qb/processedGLM/", "glm.out", full.names=T)
  load(fl[jobId])

### weather data
  load("~/Overwintering_18_19/temperatureAverage_yearFactor_GLM/weatherAve.Rdata")
  setnames(weather.ave, "V1", "sampleId")

### samps
  samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
  samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
  samps[nchar(month)==1, month:=paste("0", month, sep="")]
  samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
  samps[nchar(day)==1, day:=paste("0", day, sep="")]
  samps <- samps[set!="dgn"]
  samps[,Date:=date(paste(year, month, day, sep="-"))]

  samps[locality=="UA_od", locality:="UA_Ode"]

  samps <- merge(samps, weather.ave[,c("aveTemp", "sampleId")], by="sampleId")

  targetLocales <- list("FI_Aka"=list(locality="FI_Aka",   minYear=2000),
                        "DE_Bro"=list(locality="DE_Bro",   minYear=2000),
                        "VA_ch"=list(locality="VA_ch",     minYear=2014),
                        "PA_li"=list(locality="PA_li",     minYear=2000),
                        "DE_Mun"=list(locality="DE_Mun",   minYear=2000),
                        "UA_Ode"=list(locality=c("UA_Ode"), minYear=2010),
                        "TR_Yes"=list(locality="TR_Yes",   minYear=2000))

### open GDS
 genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)


###############
### windows ###
###############

### define windows
  win.bp <- 1e5
  step.bp <- 5e4
  setkey(glm.out, "chr")
  wins <- foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind", .errorhandling="remove")%dopar%{
      tmp <- glm.out[J(chr.i)]
      data.table(chr=chr.i,
                  start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                  end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
  }
  wins[,i:=1:dim(wins)[1]]
  dim(wins)

### run windows
  setkey(glm.out, chr, pos)
  psi <- glm.out[1]$locality

  win.out <- foreach(win.i=c(1:dim(wins)[1]), .errorhandling="remove")%dopar%{
    # win.i <- 44
    message(paste(win.i, dim(wins)[1], sep=" / "))
    win.tmp <- glm.out[J(data.table(chr=wins[win.i]$chr, pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]
    #win.tmp <- glm.out[J(data.table(index=wins[win.i]$start:wins[win.i]$end, key="index")), nomatch=0]
    win.tmp[,Z:=qnorm(p.lrt, 0, 1)]
    win.tmp[,rnpZ:=qnorm(rnp.clean, 0, 1)]

    seqSetFilter(genofile, variant.id=unique(win.tmp$id),
              sample.id=samps[
                              locality==targetLocales[which(psi==names(targetLocales))][[1]]$locality
                              ][
                              year>=targetLocales[which(psi==names(targetLocales))][[1]]$minYear]$sampleId)

    af <- seqGetData(genofile, "annotation/format/FREQ")
    f.hat <- data.table(fhat=colMeans(af[[2]], na.rm=T), id=seqGetData(genofile, "variant.id"))

    win.tmp <- merge(win.tmp, f.hat, by="id")
    win.tmp[,het:=2*fhat*(1-fhat)]

    thrs <- expand.grid(sapply(c(1:9), function(x) x*10^c(-5:-1)))[,1]
    #thrs <- c(0.01)
    tmpo <- foreach(pr.i=thrs, .errorhandling="remove", .combine="rbind")%do%{
      win.tmp[!is.na(rnp.clean),
                list(pr=pr.i,
                    rnp.pr=c(mean(rnp.clean<=pr.i)),
                    rnp.binom.p=c(binom.test(sum(rnp.clean<=pr.i), length(rnp.clean), pr.i)$p.value),
                    rnChi.pr=c(mean(rnChi.clean<=pr.i)),
                    rnChi.binom.p=c(binom.test(sum(rnChi.clean<=pr.i), length(rnChi.clean), pr.i)$p.value),
                    wZa=sum(het*Z)/(sqrt(sum(het^2))),
                    wZa.p=pnorm(sum(het*Z)/(sqrt(sum(het^2))), lower.tail=T),
                    rnp.wZa=sum(het*rnpZ)/(sqrt(sum(het^2))),
                    rnp.wZa.p=pnorm(sum(het*rnpZ)/(sqrt(sum(het^2))), lower.tail=T),
                    psi=psi,
                    mean.chisq=c(mean(chisq)),
                    median.chisq=c(median(chisq)),
                    min.p.lrt=min(p.lrt),
                    min.rank=min(rp.clean),
                    min.rnp=min(rnp.clean),
                    gr=c("obs"),
                    nSNPs=.N, i=win.i,
                    start.bp=min(win.tmp$pos), stop.bp=max(win.tmp$pos), win.i=win.i,
                    invName=paste(unique(invName), collapse=";")),
                list(mod, chr, delta, locality, perm)]
    }
    tmpo
  }
  win.out <- rbindlist(win.out)
  win.out <- merge(win.out, wins, by="i", all.y=T)
  win.out[mod=="aveTemp+year_factor"][order(Z=wZa)]
  win.out[mod=="aveTemp+year_factor"][order(rnp.binom.p)]

### save
  save(win.out, file=paste("/scratch/aob2x/dest_glm_final_nested_qb/windowAnalysis/WZA_window.dt.",
                          glm.out[1]$locality, "_", glm.out[1]$perm, ".Rdata", sep=""))
