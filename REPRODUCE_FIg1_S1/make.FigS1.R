####
rm(list = ls())

#load packages
library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(gtools)
library(patchwork)
library(ggalt)
library(reshape2)
library(ggrepel)
library(sp)
library(permute)
library(RColorBrewer)
library(data.table)
library(foreach)
library(scales)

load("S1.PCA_allpops_dims12.Rdata")

chr.pcas.outer %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    fill=year,
    #label=year,
    shape =season
  )) +
  #geom_text_repel(size = 1.5, max.overlaps = 10)+
  geom_point(size = 2) +
  scale_shape_manual(values = c(21,22,23,24)) +
  scale_fill_gradientn(
    colors=c("springgreen","cyan","blue","gold","red"),
    values=rescale(c(2011,2013,2015,2016,2018))
  ) +
  theme_bw() +
  facet_grid(city~chr) ->
  multipop_plot.chr
ggsave(multipop_plot.chr, file = "multipop_plot.chr.pdf", h = 8, w = 6.5)
