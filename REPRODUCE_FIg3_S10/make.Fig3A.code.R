#### Reproduce Figure 3
####  
####  
####  
####  Part 1: Reproduce panel A

rm(list = ls())

### libraries
library(data.table)
library(foreach)
library(doMC)
registerDoMC(4)
library(SeqArray)
library(lubridate)
library(tidyverse)
library(magrittr)
library(ggVennDiagram)
library(patchwork)

load("./dat.for.panel3A.Rdata")

ggplot() +
  geom_rect(data = final.windows.pos,
            aes(xmin=start/1e6, xmax = end/1e6,
                ymin = -110, ymax = 110), 
            alpha = 0.7, fill = "gold") +
  geom_line(
    data=dat.for.plot,
    aes(
      x=pos_mean/1e6,
      y=uci.mod,
      color = paste(perm_type, metric)
    ), alpha = 0.5
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(-115,110) +
  scale_color_manual(values = c("grey50", "grey50", "blue", "purple")) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Genomic Position (Mb)") +
  ylab("Transformed P-value") +
  xlim(0,20.5) +
  geom_vline(xintercept = 2225744/1e6) +
  geom_vline(xintercept = 13154180/1e6) ->
  panel.for.fig3

ggsave(panel.for.fig3, file = "panel.for.fig3.pdf", w = 7, h = 3)
