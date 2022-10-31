### libraries
library(data.table)
library(foreach)
library(doMC)
registerDoMC(2)


### load glm.out
fl <- list.files("/project/berglandlab/thermal_glm_dest/dest_glm_final_nested_qb/windowAnalysis", "WZA_window.dt", full.names=T)

win.out.final <- foreach(fl.i=fl, .combine="rbind")%dopar%{
  #fl.i <- fl[1]
  message(fl.i)
  load(fl.i)
  win.out
}


### save
save(win.out.final, 
     file="/scratch/yey2sn/Overwintering_ms/4.GML_plots/final.window_analysis_output.nested.qb.Rdata")

### some basic summaries,
####win.out.ag <- win.out[mod=="aveTemp+year_factor",
####                      list(min.rbinom.p=min(rbinom.p),
####                           min.wZa.p=min(wZa.p)),
####                      list(locality, perm, chr=chr.x)]
####
####win.out.ag.ag <- win.out.ag[,list(pr=mean(min.wZa.p[perm==0] > min.wZa.p[perm!=0], na.rm=T)), list####(locality, chr)]
####win.out.ag.ag
####
#load(fl[jobId])
