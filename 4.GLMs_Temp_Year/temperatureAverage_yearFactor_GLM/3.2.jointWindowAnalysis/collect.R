# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

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
  fl <- list.files("/scratch/aob2x/dest_glm_morePerms_nested_qb/join_windowAnalysis/", "joint_window.dt", full.names=T)

  win.out <- foreach(fl.i=fl, .combine="rbind")%dopar%{
    #fl.i <- fl[1]
    message(fl.i)
    load(fl.i)
    win.out[thr==0.01, pa:=p.adjust(fet.p)]
    win.out[thr==0.05, pa:=p.adjust(fet.p)]

    win.out
  }




  save(win.out, win.out.ag, win.out.ag.ag, file="/project/berglandlab/alan/joint_window.Rdata")



### load and plot
  scp aob2x@rivanna.hpc.virginia.edu:~/joint_window.Rdata ~/.

  library(data.table)
  library(ggplot2)
  load("~/joint_window.Rdata")

  inv.dt <- fread("Overwintering_18_19/Inversions/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr")



  win.out.ag <- win.out[,list(n=sum(pa<.1), nl=sum(fet.p<.0005), pops=paste(unique(locality.y[fet.p<.0005]), collapse="_")), list(thr, win.i, perm, start.bp, stop.bp, chr=chr.x)]
  win.out.ag.ag <- win.out.ag[,list(n=mean(nl), lci=quantile(nl, .025), uci=quantile(nl, .975), pops=pops[1]), list(thr, win.i, perm=I(perm!=0), start.bp, stop.bp, chr)]

  g1 <-
  ggplot() +
  geom_vline(data=inv.dt, aes(xintercept=start, linetype=invName)) +
  geom_vline(data=inv.dt, aes(xintercept=stop, linetype=invName)) +
  geom_ribbon(data=win.out.ag.ag[perm==T], aes(x=start.bp/2 + stop.bp/2, ymin=lci, ymax=uci), fill="black", alpha=.5) +
  geom_point(data=win.out.ag.ag[perm==F], aes(x=start.bp/2 + stop.bp/2, y=n, color=pops)) +
  facet_grid(~chr)

  ggsave(g1, file="~/multipop_overlap.pdf", width=10, h=4)
  
