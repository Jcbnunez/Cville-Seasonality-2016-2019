
library(tidyverse)
library(reshape2)
library(magrittr)
library(fastglm)
library(vroom)

##setwd("../18.Rudman_reanalysis/")

temp.pa.20014 = vroom("/project/berglandlab/DEST_adjacent_Rudman2022_data/orchard_2014.sig_sites.csv", delim = ",", comment = "#")
all.dat = load("/project/berglandlab/DEST_adjacent_Rudman2022_data/rudman.raw.afs.Rdata")
  
rudman.afs %>%
  filter(chr == "2L" & pos == 2209048 ) %>%
  ggplot(aes(
    x= as.numeric(Timepoint),
    y=af,
    color = Cage
  )) +
  geom_line() +
  geom_point() ->
  tsp.rudman
ggsave(tsp.rudman, file = "tsp.rudman.pdf")

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

  
  
