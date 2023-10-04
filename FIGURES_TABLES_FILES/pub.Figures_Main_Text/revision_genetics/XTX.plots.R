#### XTX megaplot script



#### XTX megaplot
xtx.enrich.plot <- ggplot(data=xtx.enrich, aes(x=as.factor(thr), y=pr/(1-thr), group=chr_class)) +
  geom_line() + facet_grid(inv~chr) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1), legend.position="bottom")


xtx.overlap.ag <- xtx.overlap[,list(pr=mean(or[perm==0] > or[perm!=0])), list(chr, inv)]
xtx.overlap.plot <- ggplot(data=xtx.overlap[perm!=0], aes((or))) + 
  geom_histogram() + 
  geom_vline(data=xtx.overlap[perm==0], aes(xintercept=(or)), color="red") + 
  facet_grid(inv~chr) +
  ylab("GLM/XtXst\noverlap") + xlab("Odds ratio") +
  geom_text(data=xtx.overlap.ag, aes(x=5, y=15, label=pr)) + theme_bw()

bf.enrich.plot <- ggplot(data=bf.enrich, aes(x=as.factor(bf_thr), y=log10(pr), color=pod_class)) + geom_point() + 
  facet_wrap(inv~chr, scales="free_y") +
  ylab("BF\nenrichment") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


bf.overlap.ag <- bf.overlap[bf_thr==10,list(pr=mean(or[perm==0] > or[perm!=0])), list(chr, inv)]
bf.overlap.plot <- ggplot(data=bf.overlap[perm!=0][bf_thr==10], aes((or))) + 
  geom_histogram() + 
  geom_vline(data=bf.overlap[perm==0][class=="real"][bf_thr==10], aes(xintercept=(or)), color="red") + 
  facet_grid(inv~chr) +
  ylab("GLM/BF\noverlap") + xlab("Odds ratio") +
  geom_text(data=bf.overlap.ag, aes(x=40, y=50, label=pr)) + theme_bw()



layout<-"
  AB
  AB
  CD
  CD"

xtx.enrich.plot + xtx.overlap.plot + bf.enrich.plot + bf.overlap.plot + plot_layout(design=layout) + plot_annotation(tag_level="A")