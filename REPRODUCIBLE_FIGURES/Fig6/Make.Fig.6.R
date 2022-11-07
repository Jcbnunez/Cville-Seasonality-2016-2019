### Figure 6
### 
rm(list = ls())


library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(foreach)
library(doParallel)
library(SNPRelate)
library(SeqArray)
library(vroom)
library(patchwork)
library(fastglm)
library("treeio")
library("ggtree")
library("TDbook")
library("adegenet")
library("FactoMineR")
library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(vcfR)
library(scatterpie)

####
gff.msp.300 = vroom("gff.msp.300.dat.txt", col_names = F, delim = "\t")
for(i in 1:dim(gff.msp.300)[1] ){
  tmp <- gff.msp.300[i,]
  
  tmp.str = str_split(tmp, pattern = "\\|")
  
  isoform = gsub('.+(isoform .).+','\\1',tmp.str[[9]][3])
  
  gff.msp.300$isoform[i] = isoform
  
}
gff.msp.300$isoform %>% unique()

gff.msp.300 %<>% 
  mutate(rank = case_when(isoform == "isoform B" ~ -0.001,
                          isoform == "isoform D" ~ -0.002,
                          isoform == "isoform E" ~ -0.003,
                          isoform == "isoform F" ~ -0.004,
                          isoform == "isoform G" ~ -0.005,
                          isoform == "isoform H" ~ -0.006,
                          isoform == "isoform I" ~ -0.007,
                          isoform == "isoform J" ~ -0.008,
                          isoform == "isoform K" ~ -0.009,
                          isoform == "isoform L" ~ -0.010,
                          isoform == "isoform M" ~ -0.011
  ))


#############
load("FST.data.fig.Rdata")

fst.dat.merged.all %>%
  filter(WS == "W_5000.S_1000") %>%
  filter(CHROM == "2L") ->
  fst.all.2L

fst.all.2L %>%
    filter(BIN_START > 4.96e6 &  BIN_END < 5.35e6 ) %>%
    ggplot(aes(
      x=(BIN_START+BIN_END)/2,
      y=WEIGHTED_FST,
      color = samp,
      #linetype = kar
    )) +
    geom_vline(xintercept = 5192177, color = "red" ) +
    geom_hline(yintercept = 0.6, color = "blue", linetype = "dashed" ) +
    #scale_color_manual(values = colors) +
    ylab(expression(paste("Wier and Cockerham ", F[ST] ) )) +
    theme_classic() +
    xlab("Genomic Position (Mb)") +
    geom_line() -> fst.dat.tsp.win
  
###### Panel B
fst.all.2L %>%
  filter(BIN_START > 4.96e6 &  BIN_END < 5.35e6 ) %>%
  arrange(-WEIGHTED_FST) %>%
  filter(WEIGHTED_FST > 0.6)

ggplot() +
  geom_rect(data = gff.msp.300,
            aes(xmin=X4, xmax = X5,
                ymin = rank-0.0001, ymax = rank+0.0001), 
            alpha = 0.7, fill = "red") +
  geom_text(
    data = slice_head(group_by(gff.msp.300, isoform)),
    aes(x=5.07e6, y = rank,
        label = isoform
    ), size = 2.5 ) +
  theme_classic() +
  xlab("Genomic Position (Mb)") +
  #geom_vline(xintercept = (5170001+5180000 )/2, color = "blue" ) +
  geom_vline(xintercept = (5190001+5200000 )/2, color = "red" ) +
  xlim(4.96e6,5.35e6) ->
  msp300.arc

ggsave((fst.dat.tsp.win/msp300.arc), file = "pi.fst.dat.5.pdf", w = 3, h =4)

##### Panel C
##### 

load("msp300.target.Rdata")

message("see all msp300 isoform annotations")
all.isoforms.annot.msp300


msp300.target %>%
  filter(set %in% c("CvilleSet", "DrosEU", "DrosRTEC")) %>%
  summarise(mean.af = mean(af, na.rm = T))

msp300.target %>%
  filter(set %in% c("CvilleSet", "DrosEU", "DrosRTEC")) %>%
  group_by(continent) %>%
  summarise(mean.af = mean(af, na.rm = T))

msp300.target %>%
  filter(sampleId == "SIM") 


load("./nasa_power.weather.mod.Rdata")
names(weather.ave)[1] = "sampleId"

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

weather.ave %>%
  filter(mod == 2) %>%
  dplyr::select(sampleId, temp.max) -> tmax


msp300.target %>%
  left_join(tmax) %>%
  filter(locality == "VA_ch") %>%
  filter(year %in% 2016:2018) %>%
  ggplot(aes(
    x=temp.max,
    y=af)
  ) + 
  theme_bw() +
  geom_smooth(method = "lm", fill = "grey", alpha = 0.3, color = "black") +
  geom_point(aes(shape=as.factor(year)), size = 3) ->
  msp300.exa

ggsave(msp300.exa, file = "msp300.exa.pdf", w = 4, h =3)

####
msp300.target %>%
  left_join(tmax) %>%
  filter(locality == "VA_ch") %>%
  filter(year %in% 2016:2018) %>%
  ggplot(aes(
    x=as.Date(collectionDate, format = "%m/%d/%Y"),
    y=af,
    color = as.factor(year))
  ) + 
  theme_bw() +
  geom_point() +
  geom_smooth(method = "lm") ->
  #geom_smooth(method = "lm", fill = "grey", alpha = 0.3, color = "black") +
  #geom_point(aes(shape=as.factor(year)), size = 3) ->
  msp300.exa.time

ggsave(msp300.exa.time, file = "msp300.exa.time.pdf", w = 4, h =3)


### model
#5170001-5200000

msp300.target %>%
  left_join(tmax) %>%
  filter(locality == "VA_ch") %>%
  filter(year %in% 2016:2018) %>%
  .[complete.cases(af),] -> dat.msp300.5192177

dat <- dat.msp300.5192177
glm.method <- 0
X.year.var <- model.matrix(~as.factor(year)+temp.max, dat)

t1.year.var <- fastglm(x=X.year.var, y=dat$af_nEff,
                       family=binomial(), weights=dat$nEff, 
                       method=glm.method)
summary(t1.year.var)

####
####
X.each.year.var <- model.matrix(~temp.max, dat)

dat <- filter(dat.msp300.5192177, year == 2016)
t2016.year.var <- fastglm(x=X.each.year.var, y=dat$af_nEff,
                          family=binomial(), weights=dat$nEff, 
                          method=glm.method)
summary(t2016.year.var)
###
dat <- filter(dat.msp300.5192177, year == 2017)
X.each.year.var <- model.matrix(~temp.max, dat)
t2017.year.var <- fastglm(x=X.each.year.var, y=dat$af_nEff,
                          family=binomial(), weights=dat$nEff, 
                          method=glm.method)
summary(t2017.year.var)
###
dat <- filter(dat.msp300.5192177, year == 2018)
X.each.year.var <- model.matrix(~temp.max, dat)
t2018.year.var <- fastglm(x=X.each.year.var, y=dat$af_nEff,
                          family=binomial(), weights=dat$nEff, 
                          method=glm.method)
summary(t2018.year.var)


###### Tree plot
###### Tree plot
###### Tree plot
###### Tree plot
###### Tree plot


nwk <- "msp.tree/msp300.treefile"
tree <- read.tree(nwk)
#genotype <- as.data.frame(vroom("tree.metadat.txt"))
#genotype.tip = as.data.frame(genotype[,c("kar")])
#row.names(genotype.tip) = paste(genotype$tip, ".0", sep = "")
dna.seq <- fasta2DNAbin("Msp.300.haps.al.polyOnly.fasta", quiet=FALSE, chunkSize=10, snpOnly=FALSE)
bin.dna <- DNAbin2genind(dna.seq , pop=NULL, exp.char=c("a","t","g","c"), polyThres=1/100)

bin.dna@tab %>%
  as.data.frame %>% 
  .[,seq(from=1, to = dim(.)[2], by = 2 )] ->
  simplified.bin.DNA

ggplot(tree, aes(x, y), 
       branch.length="none"
) + 
  geom_tree() + 
  theme_tree() +
  geom_tiplab(size=2, align=TRUE, linesize=.5) + 
  geom_nodelab(aes(x=branch, label=label), vjust=-.5, size=3) -> tree.plot


msaplot(tree.plot, dna.seq, 
        window=c(139-15, 139+15),  
        color = c("#C77CFF", "#619CFF", "#F8766D", "#00BFC4"),
        offset = 5) -> dna.tree

ggsave(dna.tree, file = "dna.tree.tree.plot.pdf") 

####
#### Map

load("msp.300.wolrd.plot.dat.Rdata")

msp.300.dat %>%
  filter(set %in% c("CvilleSet", "DrosEu", "DrosRTEC", "dgn" )) %>%
  group_by(continent) %>%
  summarise(af = mean(af))

msp.300.dat %>%
  filter(set %in% c("CvilleSet", "DrosEu", "DrosRTEC", "dgn" )) %>%
  separate(locality,
           into= c("region", "cit")) %>%
  .[complete.cases(af),] %>%
  group_by(region) %>%
  summarise(lat = mean(lat), long = mean(long), af = mean(af)) %>%
  mutate(B= round(af*100), A=round((1-af)*100)) ->
  msp.300.dat.plot

ggplot(data = world) +
  geom_sf(fill= "antiquewhite", alpha = 0.8) +
  coord_sf(xlim = c(-125.15, 60.00), ylim = c(-38.00, 65.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue")) + 
  geom_scatterpie(aes(x=long, y=lat,  r=2.5), data=msp.300.dat.plot,
                  cols=LETTERS[1:2], 
                  color="black") +
  xlab("Lon") + 
  ylab("Lat") + 
  ggtitle("D. melanogaster Msp300 mutation") + 
  theme(legend.position = "none") -> mel.world

ggsave(mel.world, file ="mel.world.pdf")

##### DPGP
######
######
######
######


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

