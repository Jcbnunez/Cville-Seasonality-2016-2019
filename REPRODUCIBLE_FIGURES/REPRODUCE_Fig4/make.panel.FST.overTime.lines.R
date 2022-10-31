### FST 
### 
### 
rm(list = ls())

setwd("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/REPRODUCE_Fig4")

library(tidyverse)
library(reshape2)
library(magrittr)
library(foreach)
library(viridis)
library(data.table)
library(ggridges)
library(broom)
library(ade4)

###save seas data
load("./FST.seasonal.Rdata")

fst.dat.EC$SNP.set %>% unique() -> sets
fst.dat.EC %<>%
  mutate(delta_y = abs(year1-year2))

###
###
cities_to_choose <- c("Munich", "Akaa", "Broggingen", "Odesa", "Charlottesville")

fst.dat.EC %>%
  filter(pop1 %in% cities_to_choose) %>%
  group_by(pop1, delta_y, SNP.set) %>%
  summarise(mean.fst = mean(FST) ) %>%
  ggplot(aes(
    x=delta_y,
    y=mean.fst,
    color=SNP.set 
  )) + 
  geom_line(size = 1.3) +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~pop1, scales = "free", nrow = 1) +
  scale_color_viridis(discrete = TRUE, option = "D") 
  






