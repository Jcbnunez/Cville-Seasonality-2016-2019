library(tidyverse)
library(vroom)
library(magrittr)
library(patchwork)

fst.2l.sim <- vroom("/Users/jcbnunez/Downloads/sim.mut.haps.fst.windowed.weir.fst")

fst.2l.sim %>%
  ggplot(aes(
    x=(BIN_START+BIN_END)/2,
    y=MEAN_FST
  )) + 
  geom_vline(xintercept = 4996892, color = "red") +
  geom_line() + 
  ggtitle("Simulans FST - T/T vs G/G homozygous") ->
  all.2l

fst.2l.sim %>%
  ggplot(aes(
    x=(BIN_START+BIN_END)/2,
    y=MEAN_FST
  )) + 
  geom_line() +
  xlim(4996892-1e5,4996892+1e5) +
  geom_vline(xintercept = 4996892, color = "red") +
  ggtitle("Simulans FST - T/T vs G/G homozygous: zoom in")->
    zoom.2l
  
  (all.2l/zoom.2l)

  
  