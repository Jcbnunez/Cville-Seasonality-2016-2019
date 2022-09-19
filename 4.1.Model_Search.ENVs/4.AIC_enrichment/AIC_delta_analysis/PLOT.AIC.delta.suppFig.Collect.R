### download
#  system("scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/environmental_ombibus_global/bestAIC.global.Rdata ~/.")

## collect

rm(list = ls())

### libraries
library(data.table)
#library(gdata)
library(foreach)
library(doMC)
library(reshape2)
registerDoMC(5)
library(tidyverse)
library(magrittr)
###
###
###

files.loc <- list.files("/scratch/yey2sn/Overwintering_ms/4.2.env.omibus.mods/out.aic.folder")

collect.aic = foreach(i=1:length(files.loc), .combine = "rbind")%do%{
  
  tmp <- get(load(paste("/scratch/yey2sn/Overwintering_ms/4.2.env.omibus.mods/out.aic.folder",
             files.loc[i], sep = "/")))
  tmp
  
}


collect.aic %>%
  filter(perm == 0) %>%
  #filter(cluster == "5.Cville") %>%
  ggplot(aes(
    x= as.numeric(cor.mods)^2,
    y = log10(mean.AICd.1),
    color = rnp.trsh,
    shape = mods.rank.typ
  )) +
  geom_point() +
  geom_smooth(se = F, #method = "lm"
              ) +
  facet_grid(chr~invName.AIC+cluster) ->
  plot.aic.diff

ggsave(plot.aic.diff, file = "plot.aic.diff.pdf", h = 6, w = 9)


collect.aic %>%
  filter(perm == 0) %>%
  #filter(cluster == "5.Cville") %>%
  ggplot(aes(
    x= mods.rank.diff,
    y = log10(mean.AICd.1),
    color = rnp.trsh,
    shape = mods.rank.typ
  )) +
  geom_point() +
  geom_smooth(se = F, #method = "lm"
              ) +
  facet_grid(chr~invName.AIC+cluster) ->
  plot.aic.Rdiff

ggsave(plot.aic.Rdiff, file = "plot.aic.Rdiff.pdf", h = 6, w = 9)


collect.aic %>%
  filter(perm == 0 & mods.rank.typ == "year") %>%
  #filter(cluster == "5.Cville") %>%
  ggplot(aes(
    x= mods.rank.diff,
    y = (mean.AICd.1),
    color = rnp.trsh,
    shape = mods.rank.typ
  )) +
  geom_point() +
  geom_smooth(se = F, #method = "lm"
  ) +
  facet_grid(chr~invName.AIC+cluster) ->
  plot.aic.Rdiff.y

ggsave(plot.aic.Rdiff.y, file = "plot.aic.Rdiff.y.pdf", h = 6, w = 9)

