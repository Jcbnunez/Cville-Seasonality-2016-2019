### Collect and plot 
### 

rm(list = ls())

library(tidyverse)
library(magrittr)
library(data.table)
library(ggrepel)
library(FactoMineR)
library(factoextra)
library(corrplot)
#library(gdata)
library(foreach)
library(doMC)
registerDoMC(5)
library(patchwork)

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))


####
wd_core = "/project/berglandlab/alan/environmental_ombibus_global/"

### Part 1 -- build the super AIC plot
load( paste(wd_core, "bestAIC.global.Rdata", sep = "/")) 
### p_lrt ==> lrt between model vs year
### min vs year aic.

### summarize
o.ag.ag <- o.ag[,
                list(.N),
                list(perm, cluster, mod, var, chr, inv=invName!="none")]

o2 <- o.ag.ag[,list(N, prop=N/sum(N), totalN=sum(N), mod, var, chr, inv),
              list(perm, cluster)]

o2.ag <- o2[,list(prop.real=prop[perm==0], totalN=totalN[perm==0],
                  prop.perm.mu=median(prop[perm!=0]),
                  prop.perm.lci=quantile(prop[perm!=0], .01),
                  prop.perm.uci=quantile(prop[perm!=0], .99),
                  prop.perm.med=median(prop[perm!=0]),
                  prop.rr=median(log2(prop[perm==0]/prop[perm!=0])),
                  prop.sd=sd(log2(prop[perm==0]/prop[perm!=0]))),
            list(chr, inv=inv, mod, var, cluster)]
o2.ag[,rr:=prop.real/prop.perm.mu]
o2.ag[,en:=(prop.real-prop.perm.mu)/prop.perm.mu]


o2.ag[order(-prop.real)][]
o2.ag[prop.real>(prop.perm.uci)][order(rr)]
o2.ag[prop.rr-2*prop.sd>0]
o2.ag[,p:=pnorm(0, prop.rr, prop.sd)]

o2.ag[cluster=="1.Europe_W"][order(-rr)][1:10]
o2.ag[cluster=="2.North_America_W"][order(-rr)][1:10]
o2.ag[cluster=="2.North_America_E"][order(-rr)][1:10]
o2.ag[cluster=="3.Europe_E"][order(-rr)][1:10]
o2.ag[cluster=="5.Cville"][order(-rr)][1:10]

### a bit more manipulation
o2.ag[,sig:=log2(prop.real/prop.perm.uci)>0 | log2(prop.real/prop.perm.lci)<0]

o2.rank <- o2.ag[chr=="2L"][inv==T][cluster=="5.Cville"]
o2.rank[,modRank:=rank(rr)]
setkey(o2.rank, mod, var)
setkey(o2.ag, mod, var)
o2.ag <- merge(o2.ag, o2.rank[,c("mod", "var", "modRank")])
o2.ag[,inv:=ifelse(inv, "Inside Inversion", "Outside Inversion")]
o2.ag[,inv:=factor(inv, levels=c("Inside Inversion", "Outside Inversion"))]

save(o2.ag, file="~/o2.globalOmnibus.Rdata")

### plot
system("scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/environmental_ombibus_global/o2.globalOmnibus.Rdata ~/.")
load("~/o2.globalOmnibus.Rdata")

o2.ag[!is.na(sig)][order(!sig)] %>%
  filter(chr == "2L" & inv == "Inside Inversion", sig == TRUE) %>%
  group_by(cluster) %>% slice_max(prop.rr, n=1) %>% as.data.frame() %>%
  left_join(sets)


aic.en.plot <-
  ggplot(data=o2.ag[!is.na(sig)][order(!sig)], aes(color=sig)) +
  geom_point(aes(x=modRank, y=(rr), shape=inv), position=position_dodge2(width=.5), size=2) +
  geom_linerange(aes(x=modRank, ymin=(prop.real/prop.perm.uci), ymax=(prop.real/prop.perm.lci)),
                 position=position_dodge2(width=.5)) +
  facet_grid(cluster~chr+inv) +
  theme_bw() +
  geom_hline(aes(yintercept=1)) + xlab("Model") + theme(axis.text.x=element_blank()) +
  scale_y_continuous(trans = "log2")

layout <- "
  CC
  CC
  "

mega <-
  aic.en.plot +
  plot_layout(design=layout, guides = "collect") +
  plot_annotation(tag_level="A")

ggsave(mega, file="./mega_AIC.omnibus.global.png", height=6, width=12)



o2.ag[!is.na(sig),
      list(.N, nSig=sum(sig), nClust=length(unique(cluster[sig==T])), nMod=length(unique(mod[sig==T]))),
      list(var, mod, chr, inv)][order(nSig)][nSig>1]

o2.ag[sig==T & var=="humidity.ave"][chr=="2L"]


o2w <- dcast(data=o2.ag, mod+var+chr+inv ~ cluster, value.var=c("rr", "sig"))
setnames(o2w, names(o2w), gsub("[0-9]\\.", "", names(o2w)))

ggplot(data=o2w) +
  geom_point(aes(x=rr_Cville, y=rr_Europe_W, color=as.factor(sig_Cville & sig_Europe_W))) +
  facet_grid(inv~chr)

o2w[sig_Cville & sig_Europe_W]











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
  theme(axis.text.y=element_blank(), axis.text.x=element_text(size=16), axis.title.x=element_text(size=18)) +
  ylab("Days prior to sampling") + xlab("") +
  coord_flip()
defn




ggplot(data=o2.ag, aes(x=modRank, y=rr, group=cluster, color=cluster)) +
  geom_point() +
  geom_line() +
  facet_grid(inv~chr)





o2w <- dcast(o2.ag, mod+var+chr+inv~cluster, value.var="rr")
setnames(o2w, names(o2w), gsub("[0-9]\\.", "", names(o2w)))

ggplot(data=o2w, aes(x=North_America_E, y=Cville)) + geom_point()

fisher.test(table(o2w$North_America_E>1, o2w$Cville>1, o2w$chr)[,,1])


ggplot(data=o2.ag, aes(x=en, y=rr, group=cluster, color=cluster)) +
  geom_point() +
  facet_grid(inv~chr)