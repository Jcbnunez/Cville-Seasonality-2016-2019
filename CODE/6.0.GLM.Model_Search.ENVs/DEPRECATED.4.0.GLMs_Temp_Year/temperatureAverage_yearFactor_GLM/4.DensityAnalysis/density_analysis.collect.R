# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

#args = commandArgs(trailingOnly=TRUE)
#jobId=as.numeric(args[1])
#
## jobId=1

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)

### load glm.out
  fl <- list.files("/scratch/aob2x/summarized_dest_glm_nested/densityAnalysis/", ".Rdata", full.names=T)

  dens.ag <- foreach(fl.i=fl, .combine="rbind")%dopar%{
    #fl.i <- fl[1]
    message(fl.i)
    load(fl.i)
    dens.ag
  }

### save
  save(dens.ag, file="~/density_analysis_output.Rdata")













### some basic summaries,
  win.out.ag <- win.out[mod=="aveTemp",
                        list(min.rbinom.p=min(rbinom.p),
                              min.wZa.p=min(wZa.p)),
                        list(locality, perm, chr=chr.x)]

  win.out.ag.ag <- win.out.ag[,list(pr=mean(min.wZa.p[perm==0] > min.wZa.p[perm!=0], na.rm=T)), list(locality, chr)]
  win.out.ag.ag


  load(fl[jobId])
