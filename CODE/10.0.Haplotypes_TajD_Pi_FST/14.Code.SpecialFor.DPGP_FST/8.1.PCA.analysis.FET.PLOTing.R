#### PCA analysis
#### 
#### 
#### 

library(vcfR)
library(adegenet)
library(vroom)
library(FactoMineR)
library(tidyverse)
library(magrittr)
library(reshape2)
library(patchwork)

dgrp.meta = vroom("dgrp.meta.txt")
names(dgrp.meta)[1] = "sampleId"
dgrp.meta %<>% mutate(set = "DGRP")
dpgp.meta = vroom("dpgp.meta.txt")
names(dpgp.meta)[1] = "sampleId"
dpgp.meta %<>% mutate(set = "DPGP")
rbind(dgrp.meta, dpgp.meta) -> joint.meta

merged.vcf="DPGP3.DGRP.2L.merged.flt.recode.vcf.gz"
#dgrp="/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz"
#dpgp="/project/berglandlab/DPGP3_chr2L/DPGP3.DGRP.2L.merged.vcf.gz"

vcf.merged<- read.vcfR(merged.vcf, verbose = FALSE)
#vcf.dpgp<- read.vcfR(dpgp, verbose = FALSE)

x.merged <- vcfR2genlight(vcf.merged)
#x.dpgp <- vcfR2genlight(vcf.dpgp)

####
x.merged.tab <- tab(x.merged)
x.merged.tab[1:10,1:10]

#x.dpgp.tab <- tab(x.dpgp)
#x.dpgp.tab[1:10,1:10]
###
dgrp.sel = joint.meta$sampleId[which(joint.meta$`In(2L)t` != "INV/ST")]
x.merged.tab[which(rownames(x.merged.tab) %in% dgrp.sel),] -> x.merged.tab.s
x.merged.tab.s[1:10,1:10]

###
data.frame(
  sampleId= row.names(x.merged.tab.s),
tsp=x.merged.tab.s[,"2L_5192177_SNP"]) %>%
  left_join(joint.meta)-> tsp.dat


PCA(x.merged.tab.s, graph = F) -> joint.pca.obj

save(joint.pca.obj, file = "joint.pca.obj.Rdata")

joint.pca.obj$ind$coord %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) %>%
  left_join(joint.meta) %>%
  left_join(tsp.dat) %>% 
  filter(tsp %in% c(0,2)) %>%
  mutate(tsp = case_when( tsp == 0 ~ "G", tsp == 2 ~ "T") ) ->
  in2lt.global.datfoPlot
  
save(in2lt.global.datfoPlot, file = "in2lt.global.datfoPlot.Rdata")
load("in2lt.global.datfoPlot.Rdata")

  ggplot() +
    geom_jitter(
    data = in2lt.global.datfoPlot,
    size = 2.8,  
    alpha = 0.5,
    aes(
    x=Dim.1,
    y=Dim.2,
    color=as.factor(tsp),
    shape = set)) +
    theme_bw() +
    theme(legend.pos = "none") -> pca.dgrp

ggsave(pca.dgrp, file = "pca.dgrp.pdf", w = 5, h = 4)

anova(lm(Dim.1 ~ set + `In(2L)t` , data = in2lt.global.datfoPlot))
anova(lm(Dim.2 ~ set + `In(2L)t` , data = in2lt.global.datfoPlot))


#in2lt.global.datfoPlot %>%
#  dplyr::select(`In(2L)t`,  set, Dim.1, Dim.2, tsp) %>%
#  melt(id = c("In(2L)t",  "set", "tsp")) %>% 
#  ggplot(aes(
#    x=paste(`In(2L)t`),
#    y=value,
#    shape = set,
#    fill =tsp
#  )) +
#  geom_jitter(alpha = 0.3, size = 3.0, position = position_dodge(width = 0.8)) + 
#  scale_shape_manual(values = 21:22) +
#  coord_flip() +
#    facet_wrap(.~variable, ncol = 1) -> pca.dgrp.box
#
#ggsave(pca.dgrp.box, file = "pca.dgrp.box.pdf", w = 4, h = 4)


#### FET tables
in2lt.global.datfoPlot %>%
  filter(set == "DGRP") %>%
  group_by(`In(2L)t`, tsp) %>%
  summarize(N = n()) %>%
  mutate(set = "DGRP")  %>%
  group_by(set) %>%
  mutate(Ntot = sum(N) ) -> DGRP.count
in2lt.global.datfoPlot %>%
  filter(set == "DPGP") %>%
  group_by(`In(2L)t`, tsp) %>%
  summarize(N = n()) %>%
  mutate(set = "DPGP") %>%
  group_by(set) %>%
  mutate(Ntot = sum(N) )-> DPGP.count


####
dat.dgrp <- data.frame(
  "STD" = c(4, 153),
  "INV" = c(15, 7),
  row.names = c("T", "Not.T"),
  stringsAsFactors = FALSE
)
dat.dgrp

chisq.test(dat.dgrp)$expected
fisher.test(dat.dgrp) -> fet.dgrp
######

dat.dpgp <- data.frame(
  "STD" = c(71, 80),
  "INV" = c(2, 42),
  row.names = c("TSP1", "TSP2"),
  stringsAsFactors = FALSE
)
dat.dpgp

chisq.test(dat.dpgp)$expected
fisher.test(dat.dpgp) -> fet.dpgp
####
rbind(DGRP.count, DPGP.count) %>%
  as.data.frame() %>%
  ggplot(aes(
    x=`In(2L)t`,
    y= N/Ntot,
    fill = tsp
  )) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~set, ncol = 1) +
  coord_flip() ->
  bar.tsp.count
ggsave(bar.tsp.count, file = "bar.tsp.count.pdf", w = 2.5, h = 4)



###
rbind(
data.frame(pop = "DGRP", lci=fet.dgrp$conf.int[1], uci=fet.dgrp$conf.int[2], est = fet.dgrp$estimate ),
data.frame(pop = "DPGP", lci=fet.dpgp$conf.int[1], uci=fet.dpgp$conf.int[2], est = fet.dpgp$estimate ) ) %>%
  as.data.frame() %>%
  ggplot(aes(
    x=pop,
    y=log2(est),
    ymin = log2(lci),
    ymax = log2(uci)
  )) +
  geom_errorbar(width = 0.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() ->
  fet_oods

ggsave(fet_oods, file = "fet_oods.pdf", w = 2.5, h = 4)

####

ggsave(pca.dgrp+bar.tsp.count+fet_oods, file = "fig6.panel.dgrp.dpgp.pdf",
       w = 8, h = 3)


####
ggplot() +
  geom_point(
    data = in2lt.global.datfoPlot,
    size = 2.4,  
    alpha = 0.5,
    aes(
      x=Dim.1,
      y=Dim.2,
      color=as.factor(`In(2L)t`),
      shape = set)) -> pca.dgrp.inv

ggsave(pca.dgrp.inv, file = "pca.dgrp.inv.pdf", w = 6, h = 4)
