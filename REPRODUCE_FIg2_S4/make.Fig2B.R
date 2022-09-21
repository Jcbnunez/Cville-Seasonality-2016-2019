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
library(data.table)
library(tidyverse)

setwd("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/REPRODUCE_FIg2_S4")
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

##### other pops






