### 

library(tidyverse)
library(magrittr)
library(lubridate)
library(tidyverse)
library(reshape2)
library(magrittr)
library(foreach)
library(viridis)
library(data.table)
library(ggridges)
library(broom)

sanches.refusta <- read.delim("~/Documents/GitHub/Cville-Seasonality-2016-2019/20.old.papers.dat/sanches_refusta_dat.txt")
load("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/17.Allele_trajectory_analysis/haplotag_snps_AFS_pol.Rdata")

haplotag_snps_AFS_pol %>%
  group_by(sampleId, collectionDate, set, year, win, yday) %>%
  filter(set %in% c("CvilleSet", "DrosEU", "DrosRTEC") ) %>%
  filter(!is.na(af_polarized)) %>%
  summarise(Mean_haplotag = mean(af, na.rm = T),
            #sd_haplotag = sd(af, na.rm = T),
            #ci_l = ci(af, na.rm = T)[2],
            #ci_h = ci(af, na.rm = T)[3],
  ) %>%
  filter(set == "CvilleSet" & win == "inv") %>%
  filter(!is.na(win)) %>%
  as.data.frame() ->
  Cville_haplotags_for_viz

###get month
Cville_haplotags_for_viz %>%
  group_by(sampleId) %>%
  mutate(collectionDate = as.Date(collectionDate, format = c("%Y-%m-%d"))) %>%
  mutate(month.num =month(collectionDate)  ) %>%
  dplyr::select(Year=year, mont.num = month.num, in2lt = Mean_haplotag) %>%
  mutate(paper = "this_study") ->
   this.study

sanches.refusta  %>%
  group_by(Year, Month) %>%
  mutate(mont.num = which(month.abb == Month)) %>%
  dplyr::select(Year=Year, mont.num = mont.num, in2lt = in2lt, sampleId = Month, paper = paper) %>%
  mutate(in2lt = in2lt/100)->
  sanches.refusta.mod

rbind(this.study, sanches.refusta.mod) -> joint.dat

joint.dat %>%
  ggplot(aes(
    x=mont.num,
    y=in2lt,
      #shape= as.factor(Year),
      shape = paper
  )) +
  geom_point(aes(fill = paper, group = paper), color = "black", size = 2) +
  geom_smooth(#color = "black",
              aes(color = paper, group = paper)
              ) +
  scale_shape_manual(values = c(21, 22)) +
  theme_bw() +
  theme(legend.position = "bottom")

###
###
joint.dat %>%
  group_by(mont.num, paper) %>%
  summarise(mean.in2lt = mean(in2lt)) %>%
  reshape2::dcast(mont.num~paper) %>%
  .[complete.cases(.),] %>%
  cor.test(~Sanches_Refusa_1990 + this_study, data = .)

### Linear models

joint.dat %>% group_by(paper) %>%
  do(mod.fit = lm(in2lt~as.numeric(mont.num), data = .)) %>%
  tidy(., mod.fit)


