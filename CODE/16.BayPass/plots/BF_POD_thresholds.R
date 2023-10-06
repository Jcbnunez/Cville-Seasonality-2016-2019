# ijob -A berglandlab -c10 -p standard --mem=80G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  .libPaths(c("/scratch/aob2x/biol4559/packages_temp/", .libPaths())); .libPaths()

  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(10)
  library(ggplot2)

### get file list
  fn <- list.files(path="/standard/vol186/bergland-lab/Gio", pattern="tempmaxpod")
  fn <- fn[grepl("betai", fn)]

  setwd("/standard/vol186/bergland-lab/Gio")
  bf.sim <- foreach(fn.i=fn)%dopar%{
    # fn.i <- fn[1]
    message(fn.i)
    tmp <- fread(fn.i)
    tmp[,pod:=as.numeric(tstrsplit(fn.i, "_")[[2]])]
    tmp[,rep:=as.numeric(tstrsplit(fn.i, "_")[[3]])]
    return(tmp)
  }
  bf.sim <- rbindlist(bf.sim)
  setnames(bf.sim, "BF(dB)", "bf_db")
  bf.sim.ag <- bf.sim[,list(bf_db.mean=mean(bf_db), bf_db.median=median(bf_db)), list(MRK, pod)]


  bf.sim.thr <- bf.sim.ag[,list(bf_db.mean=quantile(bf_db.mean, c(.95, .99, .999)),
                                  bf_db.median=quantile(bf_db.median, c(.95, .99, .999)),
                                  thr=c(.95, .99, .999)), list(pod)]
  bf.sim.thr

### save
  save(bf.sim.thr, file="~/bayPass_output/bf_threshold.Rdata")
