
###
scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/environmental_ombibus_global_permV2/glmEnrich/temp.propMax\;1\;3.Europe_E.glm_enrich.Rdata ~/glmEnrich.Rdata


### library
  library(ggplot2)
  library(data.table)

### data
#  setwd("/project/berglandlab/alan/environmental_ombibus_global_permV2/glmBins")
  load("~/glmEnrich.Rdata")

### plot
  glm.bins[,list(n=sum(N)), list(perm, perm_method)]
  ggplot(data=glm.bins[perm<10], aes(x=-log10(p_lrt.low), y=N, group=perm, color=as.factor(perm!=0))) + geom_line() + facet_grid(perm_method~chr+inv)
