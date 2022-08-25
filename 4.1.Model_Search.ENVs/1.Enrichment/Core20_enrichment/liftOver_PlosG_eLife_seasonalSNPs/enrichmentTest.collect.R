# ijob -A berglandlab_standard -c5 -p standard --mem=20G
# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)

### load data
  fl <- list.files("/scratch/aob2x/temp_seasonal_SNP_overlap/", ".Rdata", full.names=T)

  o <- foreach(fl.i=fl)%do%{
    load(fl.i)
    o
  }
  o <- rbindlist(o)

  o[,perm:=tstrsplit(psi, "_")[[3]]]
  o[,locality:=paste(tstrsplit(psi, "_")[[1]], tstrsplit(psi, "_")[[2]], sep="_")]


  o.ag <- o[,list(or=(TT[perm==0]/TF[perm==0])/(TT[perm!=0]/TF[perm!=0]), perm=perm[perm!=0]), list(paper, i, chr, locality)]
  o.ag[locality=="VA_ch"][chr=="2L"][i==0.01]

  o.ag.ag <- o.ag[,list(or=mean(log2(or)), sd=sd(log2(or))), list(paper, i, chr, locality)]

  save(o, o.ag.ag, file="~/old_papers_enrichment.Rdata")

### plot data
  load("~/old_papers_enrichment.Rdata")

  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(viridis)

  ggplot(data=o[chr=="genome"], aes(x=log10(i), y=log2(or), group=perm, color=as.factor(perm==0))) +
  geom_line() +
  facet_grid(paper~locality)
