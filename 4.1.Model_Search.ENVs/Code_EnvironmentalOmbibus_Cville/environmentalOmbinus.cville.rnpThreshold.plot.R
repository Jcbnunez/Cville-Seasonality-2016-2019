# ijob -A berglandlab_standard -c4 -p dev --mem=20G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


### libraries
  library(data.table)
  #library(gdata)
  library(foreach)
  library(doMC)
  registerDoMC(5)


### guide file
  load("/project/berglandlab/alan/environmental_ombibus/mod_var.Rdata")

###
  o.rnp <- foreach(mvi=mod_var, .combine="rbind")%dopar%{
    message(mvi)
    #mvi <- mod_var[67]
    outDir <- paste("/project/berglandlab/alan/environmental_ombibus/", mvi, sep="")

    load(file=paste(outDir, "/", mvi, ".rnp_summary.Rdata", sep=""))

    o.rnp.ag[,mod_var:=mvi]
    o.rnp.ag
  }


  o.rnp.ag <- o.rnp[,list(rr.ave=mean(log2((pr[perm!=0])/(thr))),
                          rr.sd=sd(log2((pr[perm!=0])/(thr))),
                          rr=log2(pr[perm==0]/thr)),
                     list(chr, inv, thr, mod_var)]

  save(o.rnp.ag, o.rnp, file="/project/berglandlab/alan/environmental_ombibus/mov_var.rnp.ag.Rdata")



###
  system("scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/environmental_ombibus/mov_var.rnp.ag.Rdata ~/.")

  library(data.table)
  library(ggplot2)

  load("~/mov_var.rnp.ag.Rdata")

  o.rnp.ag[,mod:=tstrsplit(mod_var, "_")[[2]]]
  o.rnp.ag[,var:=tstrsplit(mod_var, "_")[[1]]]

  p1 <-
  ggplot(data=o.rnp.ag) +
  geom_line(aes(x=log10(thr), y=2^rr, group=interaction(inv, chr), color=inv, linetype=chr)) +
  geom_point(data=o.rnp.ag[(rr.ave+2*rr.sd)<rr], aes(x=log10(thr), y=2^rr), size=.5) +
  facet_grid(mod~var)


  ggsave(p1, file="~/rnp_bigplot.pdf", h=10, w=15)
