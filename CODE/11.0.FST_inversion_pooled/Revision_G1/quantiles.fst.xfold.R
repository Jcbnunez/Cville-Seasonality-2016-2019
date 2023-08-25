### Plot quantile set
### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(foreach)
library(doParallel)
library(viridis)

### load file 
load("fig2.fst.joint.qunt.data.for.plot.Rdata")

####
#### just cville
  city.select = c("Charlottesville")

  joint.qunt.data.for.plot %>%
  filter(pop %in% city.select) %>%
  group_by(type, year_diff, pop, subset, rnp.thresh) %>%
  summarize(
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