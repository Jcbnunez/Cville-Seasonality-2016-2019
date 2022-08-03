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


load("make.Fig2B.dat.Rdata")

o2.ag %<>%
  filter(cluster == "5.Cville") %>%
  group_by(chr, inv) %>%
  arrange(rr) %>%
  mutate(modRank = 1:n())

setDT(o2.ag)

aic.en.plot <-
  ggplot(data=o2.ag[!is.na(sig)][order(!sig)], aes(color=sig)) +
  geom_point(aes(x=modRank, y=(rr), shape=inv), position=position_dodge2(width=.5), size=2) +
  geom_linerange(aes(x=modRank, ymin=(prop.real/prop.perm.uci), ymax=(prop.real/prop.perm.lci)),
                 position=position_dodge2(width=.5)) +
  facet_grid(inv~chr) +
  theme_bw() +
  geom_hline(aes(yintercept=1)) + xlab("Model") + theme(axis.text.x=element_blank()) +
  scale_y_continuous(trans = "log2")

ggsave(aic.en.plot, file = "aic.en.plot.cville.pdf", w = 9, h = 4)
