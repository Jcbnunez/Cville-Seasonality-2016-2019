#### Make Figure 2
#### 

rm(list = ls())
# Load packages

#load packages
library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(gtools)
library(patchwork)
library(reshape2)
library(ggrepel)
library(sp)
library(permute)
library(RColorBrewer)
library(devtools)
library(lubridate)
library(vroom)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(foreach)
library(doParallel)
library(viridis)
library(corrplot)
#library(gdata)
library(foreach)
library(doMC)
registerDoMC(5)



load("make.Fig2.panelA.Rdata")

ggplot() + 
  geom_point(shape =21, 
             size = 3,
             alpha = 0.9,
             data = pca_table_df,
             aes(
               x=Dim.1,
               y=Dim.2,
               fill = year)
  ) +
  scale_fill_gradient2(low = "blue", high = "red", 
                       mid = "gold",
                       midpoint = 2016) +
  theme_bw() +
  ylab("PC2") +
  xlab("PC1") +
  facet_wrap(~chr) ->
  cville_chr_pca
cville_chr_pca
ggsave(cville_chr_pca,
       file = "cville_chr_pca.pdf",
       width = 5,
       height = 4)

### Quantify
anova(lm(Dim.1 ~ year, data = filter(pca_table_df, chr == "2L" )))
anova(lm(Dim.1 ~ year, data = filter(pca_table_df, chr == "2R" )))
anova(lm(Dim.1 ~ year, data = filter(pca_table_df, chr == "3L" )))
anova(lm(Dim.1 ~ year, data = filter(pca_table_df, chr == "3R" )))


anova(lm(Dim.1 ~ `In(2L)t`, data = filter(pca_table_df, chr == "2L" )))
anova(lm(Dim.1 ~ `In(2L)t`, data = filter(pca_table_df, chr == "2R" )))
anova(lm(Dim.1 ~ `In(2L)t`, data = filter(pca_table_df, chr == "3L" )))
anova(lm(Dim.1 ~ `In(2L)t`, data = filter(pca_table_df, chr == "3R" )))

### Panels B/C
###
###
###
###
###
load("./fig2.fst.joint.qunt.data.for.plot.Rdata")

####
#### just cville
city.select = c("Charlottesville")

joint.qunt.data.for.plot %>%
  filter(pop %in% city.select) %>%
  group_by(type, year_diff, pop, subset, rnp.thresh) %>%
  dplyr::summarize(
    mean.p = mean(perc.above.control),
    sd.p = sd(perc.above.control),
    median.p = quantile(perc.above.control, 0.5),
    q25.p = quantile(perc.above.control, 0.25),
    q75.p = quantile(perc.above.control, 0.75)
  ) %>%
  ggplot(aes(
    x=type,
    y=median.p,
    ymin = q25.p,
    ymax= q75.p,
    #color = type,
    color =rnp.thresh
  )) +
  geom_hline(yintercept = 0.50, linetype = "dashed") +
  geom_hline(yintercept = 0.90, linetype = "dashed") +
  geom_errorbar(width = 0.3, position = position_dodge(width = 0.5)) +
  geom_point( size = 2.1, position = position_dodge(width = 0.5)) +
  coord_flip() +
  scale_fill_viridis(discrete=TRUE) + 
  ggtitle("Cville") +
  ylab("Median of Percentiles (FST inv > FST controls) + IQRs") +
  xlab("Correlation to In(2L)t in DGRP") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7)) +
  facet_grid(subset~as.factor(year_diff), scales =  "free_y") ->
  collect.quant.sum.cville

ggsave(collect.quant.sum.cville, file = "collect.quant.sum.cville.pdf", h = 5, w = 6)

###

joint.qunt.data.for.plot %>%
  filter(comp.set %in%  c("space.eu.e", "space.eu.w", "space.NoA.E") ) %>%
  group_by(type, year_diff, pop, subset, rnp.thresh) %>%
  dplyr::summarize(
    mean.p = mean(perc.above.control),
    sd.p = sd(perc.above.control),
    median.p = quantile(perc.above.control, 0.5),
    q25.p = quantile(perc.above.control, 0.25),
    q75.p = quantile(perc.above.control, 0.75)
  ) %>%
  ggplot(aes(
    x=type,
    y=median.p,
    ymin = q25.p,
    ymax= q75.p,
    #color = type,
    color =rnp.thresh
  )) +
  geom_hline(yintercept = 0.50, linetype = "dashed") +
  geom_hline(yintercept = 0.90, linetype = "dashed") +
  geom_errorbar(width = 0.3, position = position_dodge(width = 0.5)) +
  geom_point( size = 2.1, position = position_dodge(width = 0.5)) +
  coord_flip() +
  scale_fill_viridis(discrete=TRUE) + 
  ggtitle("Space sets in general") +
  ylab("Median of Percentiles (FST inv > FST controls) + IQRs") +
  xlab("Correlation to In(2L)t in DGRP") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7)) +
  facet_grid(subset~pop, scales =  "free_y") ->
  collect.quant.sum.space

ggsave(collect.quant.sum.space, file = "collect.quant.sum.space.pdf", h = 5, w = 6)

###
city.select.all = c("Akaa",  "Broggingen", "Cross Plains", "Linvilla", "Munich","Odesa" , "Yesiloz" )

joint.qunt.data.for.plot %>%
  filter(pop %in% city.select.all) %>%
  group_by(type, year_diff, pop, subset, rnp.thresh) %>%
  dplyr::summarize(
    mean.p = mean(perc.above.control),
    sd.p = sd(perc.above.control),
    median.p = quantile(perc.above.control, 0.5),
    q25.p = quantile(perc.above.control, 0.25),
    q75.p = quantile(perc.above.control, 0.75)
  ) %>%
  ggplot(aes(
    x=type,
    y=median.p,
    ymin = q25.p,
    ymax= q75.p,
    shape = as.factor(year_diff),
    color =rnp.thresh
  )) +
  geom_hline(yintercept = 0.50, linetype = "dashed") +
  geom_hline(yintercept = 0.90, linetype = "dashed") +
  geom_errorbar(width = 0.3, position = position_dodge(width = 0.5)) +
  geom_point( size = 2.1, position = position_dodge(width = 0.5)) +
  coord_flip() +
  scale_fill_viridis(discrete=TRUE) + 
  ggtitle("Cville") +
  ylab("Median of Percentiles (FST inv > FST controls) + IQRs") +
  xlab("Correlation to In(2L)t in DGRP") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7)) +
  facet_grid(subset~pop, scales =  "free_y") ->
  collect.quant.sum.allcities

ggsave(collect.quant.sum.allcities, file = "collect.quant.sum.allcities.pdf", h = 9, w = 9)

### Panel D
### ### Panel D
### Panel D
### Panel D
### Panel D
### Panel D
### Panel D
### Panel D
### Panel D


#setwd("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/REPRODUCE_FIg2_S4")
load("make.Fig2B.dat.Rdata")

o2.ag %>%
  filter(cluster == "5.Cville") %>%
  filter(chr == "2L", inv == "Inside Inversion" ) %>%
  #group_by(chr, inv) %>%
  #arrange(prop.real/prop.perm.lci) %>%
  arrange(rr) %>%
  mutate(modRank.in2lt = 1:n()) %>%
  dplyr::select(mod, var,  modRank.in2lt) -> o2.ag.in2lt


left_join(o2.ag, o2.ag.in2lt) -> o2.ag.ranked.all
setDT(o2.ag.ranked.all)


######
foreach(clust = c("5.Cville", "2.North_America_I95", "1.Europe_W", "3.Europe_E" ) )%do%{
  
  message(clust)
  
  
  o2.ag.ranked.all %>%
    filter(cluster == clust) ->
    o2.ag.ranked
  
  #aic.en.plot <-
  ggplot() +
    geom_linerange(data = o2.ag.ranked[!is.na(sig)][order(!sig)], 
                   aes(x=modRank.in2lt,
                       ymin=(prop.real/prop.perm.uci), 
                       ymax=(prop.real/prop.perm.lci),
                   ),color = "grey", alpha = 0.5,
                   #position=position_dodge2(width=.5)
    ) +
    geom_linerange(data=
                     o2.ag.ranked[sig == TRUE][order(!sig)],
                   aes(x=modRank.in2lt,
                       ymin=(prop.real/prop.perm.uci), 
                       ymax=(prop.real/prop.perm.lci),
                       #color = sig
                   ),alpha = 0.7, size = 0.8, color = "blue"
                   #position=position_dodge2(width=.5)
    ) +
    geom_point(data = o2.ag.ranked[!is.na(sig)][order(!sig)],
               aes(x=modRank.in2lt, y=(rr), #shape=inv
               ), color = "grey", 
               #position=position_dodge2(width=.5), 
               size=2, alpha = 0.6) +
    geom_point(data = o2.ag.ranked[sig == TRUE][order(!sig)],
               aes(x=modRank.in2lt, y=(rr), #shape=inv
               ), color = "blue", 
               #position=position_dodge2(width=.5), 
               size=2, alpha = 0.6) +
    geom_linerange(data=
                     o2.ag.ranked[!is.na(sig)][order(!sig)][var == "temp.max"][mod == 2],
                   aes(x=modRank.in2lt,
                       ymin=(prop.real/prop.perm.uci), 
                       ymax=(prop.real/prop.perm.lci),
                       color = sig
                   ),alpha = 0.7, size = 1.5
                   #position=position_dodge2(width=.5)
    ) +
    geom_point(data = 
                 o2.ag.ranked[!is.na(sig)][order(!sig)][var == "temp.max"][mod == 2],
               aes(x=modRank.in2lt, y=(rr), #shape=inv
                   fill = sig),
               #position=position_dodge2(width=.5), 
               size=3, alpha = 1.0, shape = 21) +
    geom_linerange(data=
                     o2.ag.ranked[!is.na(sig)][order(!sig)][var == "null"],
                   aes(x=modRank.in2lt,
                       ymin=(prop.real/prop.perm.uci), 
                       ymax=(prop.real/prop.perm.lci),
                       color = sig
                   ),alpha = 0.7, size = 1.5
                   #position=position_dodge2(width=.5)
    ) +
    geom_point(data = 
                 o2.ag.ranked[!is.na(sig)][order(!sig)][var == "null"],
               aes(x=modRank.in2lt, y=(rr), #shape=inv
                   fill = sig),
               #position=position_dodge2(width=.5), 
               size=3, alpha = 1.0, shape = 22) +
    geom_linerange(data=
                     o2.ag.ranked[!is.na(sig)][order(!sig)][var == "pop_year"],
                   aes(x=modRank.in2lt,
                       ymin=(prop.real/prop.perm.uci), 
                       ymax=(prop.real/prop.perm.lci),
                       color = sig
                   ),alpha = 1.0, size = 1.5
                   #position=position_dodge2(width=.5)
    ) +
    geom_point(data = 
                 o2.ag.ranked[!is.na(sig)][order(!sig)][var == "pop_year"],
               aes(x=modRank.in2lt, y=(rr), #shape=inv
                   fill = sig),
               #position=position_dodge2(width=.5), 
               size=3, alpha = 1.0, shape = 24) +
    facet_grid(inv~chr) +
    theme_bw() +
    geom_hline(yintercept=1) + 
    xlab("Model") + 
    ggtitle(unique(o2.ag.ranked$cluster)) +
    theme(axis.text.x=element_blank()) +
    scale_y_continuous(trans = "log2") ->
    mod.plot
  
  ggsave(mod.plot, file = paste(unique(o2.ag.ranked$cluster), "models.pdf", sep = ""), w = 9, h = 4)
}


