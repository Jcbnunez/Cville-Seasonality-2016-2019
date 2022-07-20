### Generate best models table
library(ggplot2)
library(data.table)
library(patchwork)
library(ggrepel)
library(tidyverse)

load("/project/berglandlab/alan/environmental_ombibus_global/o2.globalOmnibus.Rdata")

sets <- data.table(mod=c(-1, 0, 1:11),
                   label=LETTERS[1:13],
                   year=c(NA, rep(1, 12)),
                   start=c(NA, NA, 0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(NA, NA, 7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

left_join(o2.ag, sets) %>% 
  mutate(model.name = paste(var, start, end, sep = "_")) ->
  o2.ag.annot

o2.ag.annot %>%
  group_by(cluster) %>%
  filter(chr == "2L") %>%
  filter(inv == "Inside Inversion") %>%
  slice_max(rr, n = 10, with_ties = F) ->
  best.mods.tab


write.table(best.mods.tab,
            file = "best.mods.tab.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
