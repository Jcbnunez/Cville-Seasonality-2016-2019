# ijob -A berglandlab_standard -c4 -p dev --mem=20G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

args = commandArgs(trailingOnly=TRUE)
jobId=as.numeric(args[1])

#jobId=66

load("/project/berglandlab/alan/environmental_ombibus/mod_var.Rdata")
mvi <- mod_var[jobId]

### libraries
  library(data.table)
  #library(gdata)
  library(foreach)
  library(doMC)
  registerDoMC(5)

  outDir <- paste("/project/berglandlab/alan/environmental_ombibus/", mvi, sep="")


  load(file=paste(outDir, "/", mvi, ".glmRNP.Rdata", sep=""))


  thrs <- expand.grid(sapply(c(1:9), function(x) x*10^c(-5:-2)))[,1]

  o.rnp.ag <- foreach(thr.i=thrs, .combine="rbind", .errorhandling="remove")%dopar%{
    message(thr.i)
    #thr.i <- thrs[20]


    o.rnp.ag <- glm.out[,list(pr=mean(rnp<thr.i), .N,
                                 thr=thr.i),
                      list(chr, inv, perm)]
    o.rnp.ag
  }

  save(o.rnp.ag, file=paste(outDir, "/", mvi, ".rnp_summary.Rdata", sep=""))
