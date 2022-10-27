# ijob -A berglandlab_standard -c5 -p standard --mem=20G
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
  library(tidyr)

### load glm.out
  fl <- list.files("/scratch/aob2x/summarized_dest_glm_nested", "glm.out", full.names=T)

  load(fl[jobId])

### density function
  dFun <- function(x, min, max) {
    dens <- density(x, from=min, to=max)
    c(paste("x_", dens$x, sep=""), paste("y_", dens$y, sep=""))
  }


### rank norm p density
  dens.ag <- glm.out[!is.na(rnp.clean),
                      list(d=dFun(qlogis(rnp.clean), min=-6, max=6)),
                      list(chr=chr, inv=invName!="none", mod, locality, perm)]

  dens.ag[,value:=tstrsplit(d, "_")[[2]]%>%as.numeric]
  dens.ag[,variable:=tstrsplit(d, "_")[[1]]]


  dens.ag.x <- dens.ag[variable=="x"][,-c("d", "variable"),with=F]
  dens.ag.y <- dens.ag[variable=="y"][,-c("d", "variable"),with=F]

  rnp.dens.ag.w <- dens.ag.x
  setnames(rnp.dens.ag.w, "value", "x")
  rnp.dens.ag.w[,y:=dens.ag.y$value]
  rnp.dens.ag.w[,stat:="rnp"]

### chisq density
  dens.ag <- glm.out[!is.na(rnp.clean),
                      list(d=dFun(chisq, min=0, max=60)),
                      list(chr=chr, inv=invName!="none", mod, locality, perm)]

  dens.ag[,value:=tstrsplit(d, "_")[[2]]%>%as.numeric]
  dens.ag[,variable:=tstrsplit(d, "_")[[1]]]


  dens.ag.x <- dens.ag[variable=="x"][,-c("d", "variable"),with=F]
  dens.ag.y <- dens.ag[variable=="y"][,-c("d", "variable"),with=F]

  chisq.dens.ag.w <- dens.ag.x
  setnames(chisq.dens.ag.w, "value", "x")
  chisq.dens.ag.w[,y:=dens.ag.y$value]
  chisq.dens.ag.w[,stat:="chisq"]

### combine
  dens.ag <- rbind(rnp.dens.ag.w, chisq.dens.ag.w)

### save
save(dens.ag, file=paste("/scratch/aob2x/summarized_dest_glm_nested/densityAnalysis/dens.dt.",
                        glm.out[1]$locality, "_", glm.out[1]$perm, ".Rdata", sep=""))
