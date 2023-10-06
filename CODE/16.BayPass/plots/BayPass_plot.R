### libraries
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(qvalue)

### load data
  setwd("~/bayPass_output")

  load("bf_overlap.Rdata")
  load("bf_real.Rdata")
  load("xtx_overlap.Rdata")
  load("xtx_real.Rdata")
  load("bf_threshold.Rdata")

### basic distribution across genome plots
  m.xtx[,sig:=(XtXst.q<.05 & XtXst.mean>median(m.xtx$XtXst.mean))]
  xtx.dist <- m.xtx[,list(nSig=sum(sig), .N), list(chr, inv)]
  xtx.dist[,xtx.pr:=nSig/N]; xtx.dist[order(xtx.pr)]

  bf.dist <- m.bf[,list(nSig=sum(bf_db.mean>mean(bf.sim.thr[thr==.999]$bf_db.mean)), .N), list(chr, inv)]
  bf.dist[,bf.pr:=nSig/N]; bf.dist[order(bf.pr)]

  setkey(xtx.dist, chr, inv)
  setkey(bf.dist, chr,  inv)
  dist.merge <-merge(xtx.dist, bf.dist)

  ggplot(data=dist.merge, aes(x=bf.pr, y=xtx.pr)) + geom_point()

###
  thr.use<-.999
  txt_size=3
  mean_plot.bf <- ggplot(data=m.bf, aes(x=-log10(rnp), y=bf_db.median)) + geom_hex() + facet_grid(inv~chr) +
                  geom_hline(yintercept=median(bf.sim.thr[thr==thr.use]$bf_db.mean)) +
                  geom_vline(xintercept=-log10(.05)) +
                  theme_bw() + xlab("-log10(Rank-Norm P)") + ylab("dB(BF)")+
                  theme(legend.position = 'bottom')

  mean_plot.xtx <- ggplot(data=m.xtx, aes(x=-log10(rnp), y=XtXst.mean)) + geom_hex() +
                   facet_grid(inv~chr) +
                   geom_hline(yintercept=min(m.xtx[sig==T]$XtXst.mean)) +
                   geom_vline(xintercept=-log10(.05)) +
                   theme_bw() + xlab("-log10(Rank-Norm P)") + ylab("XtXst")+
                   theme(legend.position = 'bottom')


  xtx.overlap.ag <- xtx.overlap[,list(pr=mean(or[perm==0] > or[perm!=0])), list(chr, inv)]
  xtx.overlap.plot <- ggplot(data=xtx.overlap[perm!=0], aes((or))) +
    geom_histogram() +
    geom_vline(data=xtx.overlap[perm==0], aes(xintercept=(or)), color="red") +
    facet_grid(inv~chr) + xlab("Odds ratio") +
    geom_text(data=xtx.overlap.ag, aes(x=4, y=17, label=pr), size=txt_size) + theme_bw() +
    theme(legend.position = 'bottom')

  bf.overlap.ag <- bf.overlap[bf_thr==thr.use,list(pr=mean(or[perm==0] > or[perm!=0])), list(chr, inv)]
  bf.overlap.plot <- ggplot(data=bf.overlap[perm!=0][bf_thr==thr.use], aes((or))) +
    geom_histogram() +
    geom_vline(data=bf.overlap[perm==0][bf_thr==thr.use], aes(xintercept=(or)), color="red") +
    facet_grid(inv~chr) + xlab("Odds ratio") +
    geom_text(data=bf.overlap.ag, aes(x=30, y=50, label=pr), size=txt_size) + theme_bw() +
    theme(legend.position = 'bottom')

  layout<-"
  AB
  AB
  AB
  CD
  CD"

  mega <- mean_plot.xtx + mean_plot.bf + xtx.overlap.plot + bf.overlap.plot + plot_layout(design=layout) + plot_annotation(tag_level="A")

  ggsave(mega, file=paste("~/bayPass_output/mega_baypass.", thr.use, ".pdf", sep=""), h=6, w=8.5)


  prop.table(table(m.xtx$XtXst.median>60))
  m.bf[,list(pr=mean(bf_db.median>=5.5)/(1-.999)), list(chr, inv)]


xtx.enrich <- foreach(thr.i=unique(xtx.sim.thr$thr), .combine="rbind")%do%{
  m.xtx[,list(pr=mean(XtXst.median>=mean(xtx.sim.thr[thr==thr.i]$xtx.median)), thr=(1-thr.i)), list(chr, inv)][order(pr)]

}
xtx.enrich[,en:=pr/thr]
ord <- paste(rep(c("2L", "2R", "3L", "3R"), each=2), rep(c(TRUE, FALSE), 4), sep="_")
xtx.enrich[,chr_inv:=factor(paste(chr, inv, sep="_"), levels=ord)]
ggplot(data=xtx.enrich, aes(color=as.factor(thr), group=as.factor(thr), y=en, x=chr_inv)) + geom_line()
