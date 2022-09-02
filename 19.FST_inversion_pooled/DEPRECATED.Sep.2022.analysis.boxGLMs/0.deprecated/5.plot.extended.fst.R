rm(list = ls())

library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(magrittr)
library(doMC)
registerDoMC(4)
library(SeqArray)
library(lubridate)
library(tidyverse)
library(gtools)
library(poolfstat)

####
glm.fst = get(load("glm.snps.extended.Year_to_year_object.Rdata"))
####
norm.fst = get(load("macthed.controls.Inv.extended.Year_to_year_object.Rdata"))

rbind(glm.fst, norm.fst) %>% 
  mutate(delta_y = abs(year1-year2)) %>%
  group_by(delta_y, SNP.set) %>%
  summarise(me.fst = median(FST)) %>%
  ggplot(aes(
    x=delta_y,
    y=me.fst,
    color = SNP.set
  )) +
  geom_point() +
  geom_line() ->
  fst.extended.control.glm

ggsave(fst.extended.control.glm, file = "fst.extended.control.glm.pdf")
