# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  library(tidyr)
fl <- list.files("/scratch/aob2x/summarized_dest_glm/", "lmer.fetInv.", full.names=T)

o <- foreach(fl.i=fl, .combine="rbind")%do%{
  #fl.i <- fl[1]
  load(fl.i)
  tmp <- tstrsplit(fl.i, "/")%>%last%>%tstrsplit(., "\\.")
  fet.invs[,perm:=tmp[[3]]]
  fet.invs
}

o[model=="aveTemp"][factor=="2Lt"][perm=="perm0"]

save(o, file="~/fet_summary_lmer.Rdata")


### plot
  load("~/fet_summary_lmer.Rdata")
  ggplot(data=o, aes(x=log10(i), y=or, group=perm, color=as.factor(perm=="perm0"))) +
  geom_line() +
  facet_grid(~factor)


  load("~/win_out_lmer.Rdata")

  inv.dt <- fread("Overwintering_18_19/Inversions/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr.x")

  ace.dt <- data.table(chr.x="3R", pos=13174333/2  + 13274333/2)


  mh.plot <- ggplot() +
  geom_vline(data=inv.dt, aes(xintercept=start, linetype=invName)) +
  geom_vline(data=inv.dt, aes(xintercept=stop, linetype=invName)) +
  geom_vline(data=ace.dt, aes(xintercept=pos), color="orange", size=2, alpha=.75) +
  geom_point(data=win.out[mod=="aveTemp"][pr==.05][nSNPs>50][order(-binom.p)],
            aes(x=start/2 +end/2 , y=-1*log10((binom.p)),
                color=rnp.pr), size=.95) +
  facet_grid(locality~chr.x, scales="free_x") +
  theme_bw() +
  scale_color_viridis(option="G", direction = -1) +
  theme(axis.text.x=element_text(angle=90))

  mh.plot

  ggsave(mh.plot, file="~/mh_lmer.pdf", h=4, w=8)
