### libraries
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(viridis)

### bring in analsysis that overlaps original Machado et al seaonsal SNPs
# scp aob2x@rivanna.hpc.virginia.edu:~/drosRtec_enrichment.Rdata ~/.
  load("~/drosRtec_enrichment.Rdata")
  o <- o[mod=="orig"][set=="aveTemp+year_factor"]
  o.drosRTEC <- o

### bring in analysis that re-implements Core20 model (with and without VA) using DEST + quasibinomial model
# scp aob2x@rivanna.hpc.virginia.edu:~/dest_core20.enrichment_signTest.Rdata ~/.
  load("~/dest_core20.enrichment_signTest.Rdata")

### clinal analysis
  # scp aob2x@rivanna.hpc.virginia.edu:~/drosRtec_enrichment_cline.Rdata ~/.
  load("~/drosRtec_enrichment_cline.Rdata")

### plots
  drosrtec <-
  ggplot() +
  geom_hline(yintercept=1) +
  geom_line(data= o.drosRTEC, aes(x=log10(thr), y=or, group=interaction(chr, inv), color=chr, linetype=inv)) +
  geom_point(data=o.drosRTEC[p<.05], aes(x=log10(thr), y=or, group=interaction(chr, inv), color=chr)) +
  geom_point(data=o.drosRTEC[p<.05/8], aes(x=log10(thr), y=or, group=interaction(chr, inv)), color="black", size=.25) +
  facet_grid(set~mod) +
  theme_bw() +
  theme(text=element_text(size=18)) +
  ylim(-.05, 3)

  drosrtec.cline <-
  ggplot() +
  geom_hline(yintercept=1) +
  geom_line(data= o.cline, aes(x=log10(thr), y=or, group=interaction(chr, inv), color=chr, linetype=inv)) +
  geom_point(data=o.cline[p<.05], aes(x=log10(thr), y=or, group=interaction(chr, inv), color=chr)) +
  geom_point(data=o.cline[p<.05/8], aes(x=log10(thr), y=or, group=interaction(chr, inv)), color="black", size=.25) +
  facet_grid(set~mod) +
  theme_bw() +
  theme(text=element_text(size=18)) +
  ylim(-.05, 3)





  core20 <-
  ggplot() +
  geom_hline(yintercept=1) +
  geom_line(data= o.dest_core20, aes(x=log10(thr), y=or, group=interaction(chr, inv), color=chr, linetype=inv)) +
  geom_point(data=o.dest_core20[p<.05], aes(x=log10(thr), y=or, group=interaction(chr, inv), color=chr)) +
  geom_point(data=o.dest_core20[p<.05/8], aes(x=log10(thr), y=or, group=interaction(chr, inv)), color="black", size=.25) +
  facet_grid(set~mod) +
  theme_bw() +
  theme(text=element_text(size=18)) +
  ylim(-.05, 3)


  core20.st <-
  ggplot() +
  geom_hline(yintercept=.5) +
  geom_line(data= o.dest_core20, aes(x=log10(thr), y=st.pr, group=interaction(chr, inv), color=chr, linetype=inv)) +
  geom_point(data=o.dest_core20[st.p<.05], aes(x=log10(thr), y=st.pr, group=interaction(chr, inv), color=chr)) +
  geom_point(data=o.dest_core20[st.p<.05/8], aes(x=log10(thr), y=st.pr, group=interaction(chr, inv)), color="black", size=.25) +
  facet_grid(set~mod) +
  theme_bw() +
  theme(text=element_text(size=18)) +
  ylim(0,1)




  layout <- "
  ABB
  #CC
  "
  mega <- drosrtec + core20 + core20.st + plot_layout(design=layout, guides="collect")

  ggsave(mega, file="~/enrichment_test.pdf", h=10, w=12)
