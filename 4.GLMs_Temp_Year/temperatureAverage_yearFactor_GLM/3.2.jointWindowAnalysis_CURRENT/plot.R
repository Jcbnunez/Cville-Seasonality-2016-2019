

### load and plot
  #scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/joint_window.Rdata ~/.

  
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  
  load("/project/berglandlab/alan/joint_window.Rdata")

  #inv.dt <- fread("Overwintering_18_19/Inversions/InversionsMap_hglft_v6_inv_startStop.txt")
  inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions/InversionsMap_hglft_v6_inv_startStop.txt"
  ### load suppl data
  inv.dt <- fread(inversion_map)
  
  
  setnames(inv.dt, "chrom", "chr")


  win.out = win.out[chr.x == "2L"]
  
  win.out.ag <- win.out[,list(n=sum(pa<.1), nl=sum(fet.p<.0005), pops=paste(unique(locality.y[fet.p<.0005]), collapse="_")), list(thr, win.i, perm, start.bp, stop.bp, chr=chr.x)]
  win.out.ag.ag <- win.out.ag[,list(n=mean(nl), lci=quantile(nl, .025), uci=quantile(nl, .975), pops=pops[1]), list(thr, win.i, perm=I(perm!=0), start.bp, stop.bp, chr)]

  g1 <-
  ggplot() +
  geom_vline(data=inv.dt[chr == "2L"], aes(xintercept=start, linetype=invName)) +
  geom_vline(data=inv.dt[chr == "2L"], aes(xintercept=stop, linetype=invName)) +
  geom_ribbon(data=win.out.ag.ag[perm==T][thr==0.05], aes(x=start.bp/2 + stop.bp/2, ymin=lci, ymax=uci), fill="black", alpha=.5) +
  geom_point(data=win.out.ag.ag[perm==F][thr==0.05], aes(x=start.bp/2 + stop.bp/2, y=n, color=pops), size = 2, alpha = 0.5) +
  facet_grid(~chr) +
  theme_bw()

  ggsave(g1, file="./multipop_overlap.pdf", h=2.9, w=5)
