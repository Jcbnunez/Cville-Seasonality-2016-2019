
library(vroom)
library(tidyverse)
library(patchwork)
library(reshape2)

fst.dgrp1 = vroom("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/16.Haplotypes_TajD_Pi_FST/data/FSTs.2L.dgrp.cm/DGRP.FST.W_100000.S_50000.INVvsSTD.windowed.weir.fst") %>% mutate(samp = "DGRP", WS = "W_100000.S_50000")
fst.cm1 = vroom("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/16.Haplotypes_TajD_Pi_FST/data/FSTs.2L.dgrp.cm/FST.W_100000.S_50000.INVvsSTD.windowed.weir.fst") %>% mutate(samp = "CM", WS = "W_100000.S_50000")

fst.dgrp2 = vroom("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/16.Haplotypes_TajD_Pi_FST/data/FSTs.2L.dgrp.cm/DGRP.FST.W_500000.S_100000.INVvsSTD.windowed.weir.fst") %>% mutate(samp = "DGRP", WS = "W_500000.S_100000")
fst.cm2 = vroom("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/16.Haplotypes_TajD_Pi_FST/data/FSTs.2L.dgrp.cm/FST.W_500000.S_100000.INVvsSTD.windowed.weir.fst") %>% mutate(samp = "CM", WS = "W_500000.S_100000")

fst.dgrp3 = vroom("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/16.Haplotypes_TajD_Pi_FST/data/FSTs.2L.dgrp.cm/DGRP.FST.W_10000.S_5000.INVvsSTD.windowed.weir.fst") %>% mutate(samp = "DGRP", WS = "W_10000.S_5000")
fst.cm3 = vroom("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/16.Haplotypes_TajD_Pi_FST/data/FSTs.2L.dgrp.cm/FST.W_10000.S_5000.INVvsSTD.windowed.weir.fst") %>% mutate(samp = "CM", WS = "W_10000.S_5000")

fst.dgrp4 = vroom("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/16.Haplotypes_TajD_Pi_FST/data/FSTs.2L.dgrp.cm/DGRP.FST.W_5000.S_1000.INVvsSTD.windowed.weir.fst") %>% mutate(samp = "DGRP", WS = "W_5000.S_1000")
fst.cm4 = vroom("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/16.Haplotypes_TajD_Pi_FST/data/FSTs.2L.dgrp.cm/FST.W_5000.S_1000.INVvsSTD.windowed.weir.fst") %>% mutate(samp = "CM", WS = "W_5000.S_1000")


rbind(fst.dgrp1, fst.cm1, fst.dgrp2, fst.cm2, fst.dgrp3, fst.cm3, fst.dgrp4, fst.cm4) %>%
  mutate(mid = (BIN_START+BIN_END)/2) %>%
  ggplot(aes(
    x=mid,
    y=MEAN_FST,
    color = samp
  )) +
  #geom_hline(yintercept = 0.6 , linetype = "solid", color = "blue") +
  geom_vline(xintercept = 5192177 , linetype = "solid", color = "blue") +
  geom_vline(xintercept = 2225744, linetype = "dashed") +
  geom_vline(xintercept = 13154180, linetype = "dashed") +
  geom_line(size = 1.0, alpha = 0.9) +
  facet_grid(WS~samp) +
  ggtitle("Unweighted FST") +
  ylim(-0.1, 0.70) +
  theme_bw() -> meanfst

rbind(fst.dgrp1, fst.cm1, fst.dgrp2, fst.cm2, fst.dgrp3, fst.cm3, fst.dgrp4, fst.cm4) %>%
  mutate(mid = (BIN_START+BIN_END)/2) %>%
  ggplot(aes(
    x=mid,
    y=WEIGHTED_FST,
    color = samp
  )) +
  #geom_hline(yintercept = 0.6 , linetype = "solid", color = "blue") +
  geom_vline(xintercept = 5192177 , linetype = "solid", color = "blue") +
  geom_vline(xintercept = 2225744, linetype = "dashed") +
  geom_vline(xintercept = 13154180, linetype = "dashed") +
  geom_line(size = 1.0, alpha = 0.9) +
  facet_grid(WS~samp) +
  ggtitle("Weighted FST") +
  ylim(-0.1, 0.70) +
  theme_bw() -> weightedfst


ggsave(weightedfst /meanfst , file = "fsts.in2L.plot.pdf", w = 9, h = 7)
####
##


#### APanel Figure
rbind(fst.dgrp1, fst.cm1) %>% 
  mutate(mid = (BIN_START+BIN_END)/2) %>%
  melt(id= c("CHROM", "mid", "BIN_START", "BIN_END", "N_VARIANTS", "samp",  "WS")) %>%
  filter(variable == "WEIGHTED_FST") %>%
  ggplot(aes(
    x=mid,
    y=value,
    color = samp,
    linetype = variable
  )) +
  #geom_vline(xintercept = 5192177 , linetype = "solid", color = "blue") +
  geom_vline(xintercept = 2225744, linetype = "dashed") +
  geom_vline(xintercept = 13154180, linetype = "dashed") +
  geom_line(size = 1.0, alpha = 0.9) +
  #facet_grid(WS~samp) +
  xlim(0, 21e6) +
  ggtitle("Weighted FST") +
  theme_bw()
  


  