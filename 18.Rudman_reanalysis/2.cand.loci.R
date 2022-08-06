
library(tidyverse)
library(reshape2)
library(magrittr)
library(fastglm)
library(vroom)

setwd("../18.Rudman_reanalysis/")

temp.pa.20014 = vroom("./orchard_2014.sig_sites.csv", delim = ",", comment = "#")

final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )

#####
temp.pa.20014 %>%
  filter(
    chrom == "2L",
  ) ->filt.2L.rud
  
ggplot() +
  geom_rect(data = final.windows.pos,
            aes(xmin=start, xmax = end,
                ymin = -2, ymax = 2), 
            alpha = 0.7, fill = "gold") +
  geom_point(data = filt.2L.rud, aes(
    x=pos,
    y=coef.div , 
    color = -log10(p.div)
  )) +
  facet_grid(comparison~.) +
  ylab("GLM Coefficient") +
  xlab("Genomic Position (Mb)") +
  ggtitle("Regions of Interest", subtitle = "Rudman et. al, 2022 Analysis") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept  = 5192177, color = "red") +
  geom_vline(xintercept = 2225744) +
  geom_vline(xintercept = 13154180) ->
  manhattan.plot.rudman
ggsave(manhattan.plot.rudman, file = "manhattan.plot.rudman.pdf")

  
  
