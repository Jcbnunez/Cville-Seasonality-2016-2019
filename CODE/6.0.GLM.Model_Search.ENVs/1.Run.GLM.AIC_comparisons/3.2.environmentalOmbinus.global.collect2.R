# ijob -A berglandlab_standard -c20 -p standard --mem=80G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

args = commandArgs(trailingOnly=TRUE)
jobId=as.numeric(args[1])

#jobId=100


### libraries
  library(data.table)
  #library(gdata)
  library(foreach)
  library(doMC)
  registerDoMC(16)

### load guide file
  #load(fl[1])
  #mod_var <- as.data.table(expand.grid(unique(glm.out$variable), unique(glm.out$mod), unique(glm.out$cluster)))
  #setnames(mod_var, names(mod_var), c("variable", "mod", "cluster"))
  #mod_var <- mod_var[!(variable=="null" & mod>-1)]
  #mod_var <- mod_var[!(variable=="pop_year" & mod>0)]
  #if(any(mod_var$cluster=="2.North_America_Mid")) {
  #       mod_var[cluster=="2.North_America_E", cluster:="2.North_America_I95"]
  #}
#
#
  #save(mod_var, file="/project/berglandlab/alan/environmental_ombibus_global/mod_var.redoNAM.Rdata")

  load("/project/berglandlab/alan/environmental_ombibus_global/mod_var.redoNAM.Rdata")
  mainDir <- "/scratch/aob2x/environmental_ombibus_global"
  mvi <- paste(mod_var[jobId]$variable, mod_var[jobId]$mod, mod_var[jobId]$cluster, sep=";")

### cycle through raw data and save mvi
### combine with snp identity
  load("~/Overwintering_18_19/SNP_filtering/snp_dt_25percMissing.Rdata")
  setkey(snp.dt, id)

  use <- apply(snp.dt[,c("VA_ch"), with=F],
                1, any)
  snp.dt <- snp.dt[use]
  setkey(snp.dt, id)


### split output and get summaries
  fl <- list.files("/scratch/aob2x/environmental_ombibus_global",full.names=T)
  fl <- fl[grepl("redoNAM", fl)]

  outDir <- paste("/project/berglandlab/alan/environmental_ombibus_global/", mvi, sep="")
  dir.create(file.path(outDir), showWarnings = FALSE)

  glm.out <- foreach(fl.i=fl, .errorhandling="remove")%dopar%{
    # fl.i <- fl[1]
    message(fl.i)
    load(fl.i)
    glm.out[cluster=="2.North_America_E", cluster:="2.North_America_I95"]
    glm.out <- merge(glm.out, snp.dt, by.x="variant.id", by.y="id", all.x=T)
    setkey(glm.out, variable, mod, cluster)

    tmp <- glm.out[J(mod_var[jobId])]
    tmp
  }
  glm.out <- rbindlist(glm.out)


  glm.out.ag <- glm.out[,list(rnp=rank(p_lrt)/length(p_lrt), chr, pos, inv=invName!="none"), list(perm)]

  setkey(glm.out, chr, pos, perm)
  setkey(glm.out.ag, chr, pos, perm)

  glm.out <- merge(glm.out, glm.out.ag)

  save(glm.out, file=paste(outDir, "/", mvi, ".glmRNP.Rdata", sep=""))


#  thrs <- expand.grid(sapply(c(1:9), function(x) x*10^c(-5:-2)))[,1]
#
#  o.rnp.ag <- foreach(thr.i=thrs, .combine="rbind", .errorhandling="remove")%dopar%{
#    message(thr.i)
#    #thr.i <- thrs[20]
#    o.rnp.ag <- glm.out[,list(pr=mean(rnp<thr.i), .N,
#                                 thr=thr.i),
#                      list(chr, inv, perm)]
#    o.rnp.ag
#  }
#
#  save(o.rnp.ag, file=paste(outDir, "/", mvi, ".rnp_summary.Rdata", sep=""))
#
