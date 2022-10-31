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

files.loc <- list.files("./out.aic.folder")

collect.aic = foreach(i=1:length(files.loc), .combine = "rbind")%do%{
  
  tmp <- get(load(paste("./out.aic.folder",
             files.loc[i], sep = "/")))
  tmp
  
}


collect.aic %>%
  filter(perm == 0) %>%
  ggplot(aes(
    x= as.numeric(cor.mods)^2,
    y = (mean.AICd.1),
    color = rnp.trsh,
    shape = mods.rank.typ
  )) +
  #geom_point() +
  geom_smooth(se = F, #method = "lm"
              ) +
  geom_hline(yintercept = (2)) +
  scale_y_continuous(trans='log10') +
  theme_bw() +
  theme(legend.position = "top",
        axis.text = element_text(size = 7))   +
  facet_grid(chr~invName.AIC+cluster) ->
  plot.aic.diff

ggsave(plot.aic.diff, file = "plot.aic.diff.pdf", h = 6, w = 9)


####### plot cville model
collect.aic$rnp.trsh = factor(collect.aic$rnp.trsh, levels = c("noRNP", "rnp5%", "rnp1%","rnp01%"))

collect.aic %>%
  filter(perm == 0) %>%
  filter(mods.rank.typ %in% c("null", "year")  ) %>%
  filter(best.mod == "temp.max;2;5.Cville" ) %>%
  ggplot(aes(
    x= mods.rank.typ,
    y = (med.AICd.1),
    ymin =  lci.AICd.1,
    ymax = uci.AICd.1,
    color = rnp.trsh,
    fill = rnp.trsh,
  )) +
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.8),size = 1.4) +
  geom_point(position = position_dodge(width = 0.8), size = 2, shape = 21, color = "black") +
  scale_y_continuous(trans='log10') +
  theme_bw() +
  theme(legend.position = "top",
        axis.text = element_text(size = 7))   +
  geom_hline(yintercept = (2)) +
  facet_grid(chr~invName.AIC) ->
  plot.aic.Rdiff

ggsave(plot.aic.Rdiff, file = "plot.aic.Rdiff.pdf", h = 4, w = 5)


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

