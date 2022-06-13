
### plot
system("scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/environmental_ombibus/bestAIC.Rdata ~/.")

### libraries
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(ggrepel)

### load data
  load("~/bestAIC.Rdata")
  load("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/SNP_filtering/snp_dt_25percMissing.Rdata")

  setnames(snp.dt, "id", "variant.id")
  setkey(snp.dt, variant.id)
  setkey(o.ag, variant.id)
  o.ag <- merge(o.ag, snp.dt[,c("chr", "pos", "variant.id", "invName")])

### summarize
  o.ag.ag <- o.ag[,
                  list(N=length(variant.id)),
                  list(perm, mod, var, chr, inv=invName!="none")]

  o2 <- o.ag.ag[,list(N, prop=N/sum(N), totalN=sum(N), mod, var),
                      list(perm , chr, inv)]

  o2.ag <- o2[,list(prop.real=prop[perm==0], totalN=totalN[perm==0],
                    prop.perm.mu=mean(prop[perm!=0]),
                    prop.perm.lci=quantile(prop[perm!=0], .01),
                    prop.perm.uci=quantile(prop[perm!=0], .99),
                    prop.perm.med=median(prop[perm!=0]),
                    prop.rr=mean(log2(prop[perm==0]/prop[perm!=0])),
                    prop.sd=sd(log2(prop[perm==0]/prop[perm!=0]))),
                list(chr, inv, mod, var)]
  o2.ag[,rr:=prop.real/prop.perm.mu]

  o2.ag[order(-prop.real)][]
  o2.ag[prop.real>(prop.perm.uci)][order(rr)]
  o2.ag[prop.rr-2*prop.sd>0]
  o2.ag[,p:=pnorm(0, prop.rr, prop.sd)]


### a bit more manipulation
  o2.ag[,sig:=log2(prop.real/prop.perm.uci)>0 | log2(prop.real/prop.perm.lci)<0]

  o2.rank <- o2.ag[chr=="2L"][inv==T]
  o2.rank[,modRank:=rank(rr)]
  setkey(o2.rank, mod, var)
  setkey(o2.ag, mod, var)
  o2.ag <- merge(o2.ag, o2.rank[,c("mod", "var", "modRank")])
  o2.ag[,inv:=ifelse(inv, "Inside Inversion", "Outside Inversion")]
  o2.ag[,inv:=factor(inv, levels=c("Inside Inversion", "Outside Inversion"))]

### plot

  aic.en.plot <-
  ggplot(data=o2.ag[!is.na(sig)][order(!sig)], aes(color=sig)) +
  geom_point(aes(x=modRank, y=(rr), shape=inv), position=position_dodge2(width=.5), size=2) +
  geom_linerange(aes(x=modRank, ymin=(prop.real/prop.perm.uci), ymax=(prop.real/prop.perm.lci)),
                position=position_dodge2(width=.5)) +
  facet_grid(inv~chr) +
  theme_bw() +
  geom_hline(aes(yintercept=1)) + xlab("Model") + theme(axis.text.x=element_blank()) +
  scale_y_continuous(trans = "log2") +
  geom_text_repel(data=o2.ag[!is.na(sig)][sig==T],
                  aes(x=modRank, y=(rr), label=paste(var, mod, sep=" / ")),
                  xlim=c(1,50),
                  size = 3,
                  point.padding = 1,
                  min.segment.length = 0,
                  max.time = 1, max.iter = 1e5,
                  box.padding = 1)



  layout <- "
  CC
  CC
  "


  mega <-
  aic.en.plot +
  plot_layout(design=layout, guides = "collect") +
  plot_annotation(tag_level="A")

  ggsave(mega, file="~/mega_AIC.omnibus.png", height=6, width=12)


  sets <- data.table(mod=c(-1, 0, 1:11),
                    label=LETTERS[1:13],
                    year=c(NA, rep(1, 12)),
                     start=c(NA, NA, 0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                     end=	 c(NA, NA, 7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))



  defn <-
  ggplot(data=sets) +
  geom_linerange(aes(x=mod, ymin=-start, ymax=-end)) +
  geom_text(aes(x=mod, y=5, label=label)) +
  geom_point(aes(x=mod, y=13+year)) +
  geom_line(aes(x=mod,  y=13+year)) +
  theme_minimal()  +
  theme(axis.text.y=element_blank()) +
  ylab("Days prior to sampling") + xlab("") +
  coord_flip()
  defn
