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
  #fl <- list.files("/scratch/aob2x/summarized_dest_glm/windowAnalysis/", ".Rdata", full.names=T)
  #fl <- list.files("/scratch/aob2x/summarized_dest_glm_nested_/windowAnalysis/", ".Rdata", full.names=T)
  fl <- list.files("/scratch/aob2x/dest_glm_morePerms_nested_qb/windowAnalysis/", "WZA_window.dt", full.names=T)

  win.out <- foreach(fl.i=fl, .combine="rbind")%dopar%{
    #fl.i <- fl[1]
    message(fl.i)
    load(fl.i)
    win.out
  }

### save
  #save(win.out, file="~/Overwintering_18_19/temperatureAverage_yearFactor_GLM/3.windowAnalysis/window_analysis_output.Rdata")
  #save(win.out, file="/scratch/aob2x/window_analysis_output.nested.qb.Rdata")
  save(win.out, file="/scratch/aob2x/dest_glm_morePerms_nested_qb/window_analysis_output.nested.qb.Rdata")




### some basic summaries,
  win.out.ag <- win.out[mod=="aveTemp",
                        list(min.rbinom.p=min(rbinom.p),
                              min.wZa.p=min(wZa.p)),
                        list(locality, perm, chr=chr.x)]

  win.out.ag.ag <- win.out.ag[,list(pr=mean(min.wZa.p[perm==0] > min.wZa.p[perm!=0], na.rm=T)), list(locality, chr)]
  win.out.ag.ag


  load(fl[jobId])
