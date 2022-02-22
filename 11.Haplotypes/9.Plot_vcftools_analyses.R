### Collect VCFTools analyses
### 

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(foreach)
library(doParallel)

inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions/InversionsMap_hglft_v6_inv_startStop.txt"
inv.dt <- fread(inversion_map)
setnames(inv.dt, "chrom", "chr.x")

####### Load all
####### Load all
####### Load all
####### Load all
####### Load all
####### Load all
####### Load all

pi_cm <- fread("./PI_D_AllSamps/Pi.CM.W_100000.S_50000.All.windowed.pi")
pi_cm %>%
  mutate(mid_point = (BIN_START+BIN_END)/2,
         pop = "CM",
         type = "All") %>%
  dplyr::select(PI,  mid_point, pop, type) %>%
  melt(id = c("mid_point", "pop", "type")) -> pi_cm_parsed

pi_dgrp <- fread("./PI_D_AllSamps/Pi.DGRP.W_100000.S_50000.All.windowed.pi")
pi_dgrp %>%
  mutate(mid_point = (BIN_START+BIN_END)/2,
         pop = "DGRP",
         type = "All") %>%
  dplyr::select(PI,  mid_point, pop, type) %>%
  melt(id = c("mid_point", "pop", "type")) -> pi_dgrp_parsed

D_cm <- fread("./PI_D_AllSamps/D.CM.W_100000.S_50000.All.Tajima.D")
D_cm %>%
  mutate(mid_point = BIN_START,
         pop = "CM",
         type = "All") %>%
  dplyr::select(TajimaD,  mid_point, pop, type) %>%
  melt(id = c("mid_point", "pop", "type")) -> D_cm_parsed

D_dgrp <- fread("./PI_D_AllSamps/D.DGRP.W_100000.S_50000.All.Tajima.D")
D_dgrp %>%
  mutate(mid_point = BIN_START,
         pop = "DGRP",
         type = "All") %>%
  dplyr::select(TajimaD,  mid_point, pop, type) %>%
  melt(id = c("mid_point", "pop", "type")) -> D_dgrp_parsed

##
rbind(pi_cm_parsed,pi_dgrp_parsed,D_cm_parsed,D_dgrp_parsed ) -> all_pi_d_parsed

####### Load subset
####### Load subset
####### Load subset
####### Load subset
####### Load subset

files_pi = system("ls ./Pi_D_subsets | grep 'pi'", intern = T)

pi_merge <- foreach(i=1:length(files_pi), .combine = "rbind")%do%{
info <- strsplit(files_pi[i], split = "\\.")[[1]]
tmp <- fread(paste("./Pi_D_subsets/", files_pi[i], sep = ""))
tmp %>%
  mutate(pop = info[2],
         type = info[5],
         resolution = info[3]) %>%
  dplyr::select(PI,  BIN_START, pop, type, resolution) %>%
  melt(id = c("BIN_START", "pop", "type", "resolution")) -> tmp_parsed
  return(tmp_parsed)
}

files_D = system("ls ./Pi_D_subsets | grep 'Taj'", intern = T)

D_merge <- foreach(i=1:length(files_D), .combine = "rbind")%do%{
  info <- strsplit(files_D[i], split = "\\.")[[1]]
  tmp <- fread(paste("./Pi_D_subsets/", files_D[i], sep = ""))
  tmp %>%
    mutate(pop = info[2],
           type = info[5],
           resolution = info[3]) %>%
    dplyr::select(TajimaD,  BIN_START, pop, type, resolution) %>%
    melt(id = c("BIN_START", "pop", "type", "resolution")) -> tmp_parsed
  return(tmp_parsed)
}


rbind(pi_merge,D_merge ) -> sub_pi_d_parsed

###### plot
###### plot
###### plot
###### plot
###### plot
###### plot


###### save
###### 
#save(inv.dt,all_pi_d_parsed,sub_pi_d_parsed, file = "pi_D_datforplot.Rdata")
load("./pi_D_datforplot.Rdata")

all_pi_d_parsed %>%
  ggplot(
    aes(
      x=mid_point,
      y=value,
      #color =karyo,
      linetype = pop)
  ) + 
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start), linetype="dashed") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop), linetype="dashed") +
  geom_line(alpha = 0.5) +
  ggtitle("All Samples") +
  theme_bw() +
  facet_wrap(~variable, ncol = 1, scales = "free") ->
  pi_d_plot_all

ggsave(pi_d_plot_all, file = "pi_d_plot_all.pdf")

sub_pi_d_parsed %>% 
  group_by(pop, type, resolution) %>%
  summarize(N = n())

sub_pi_d_parsed %>%
  filter(resolution == "W_100000") %>%
  ggplot(
    aes(
      x=BIN_START,
      y=value,
      color =type,
      linetype = pop)
  ) + 
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start), linetype="dashed") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop), linetype="dashed") +
  geom_line(alpha = 0.5) +
  ggtitle("Samples split by Karyotype") +
  theme_bw() +
  facet_wrap(~variable, ncol = 1, scales = "free") ->
  pi_d_plot_sub

ggsave(pi_d_plot_sub, file = "pi_d_plot_sub.pdf")

###### zoom Window
###### zoom Window
###### zoom Window
###### zoom Window
###### zoom Window
#####  ######
load("./pi_D_datforplot.Rdata")

sub_pi_d_parsed %>%
  filter(resolution == "W_10000") %>% 
  mutate(zoom_win = case_when(BIN_START > 5.5e6 & BIN_START < 7e6 ~ "win5",
                              BIN_START > 9.0e6 & BIN_START < 10.6e6 ~ "win9",
                              )) %>% 
  filter(zoom_win %in% c("win5", "win9")) %>%
  ggplot(
    aes(
      x=BIN_START,
      y=value,
      color =type,
      linetype = pop)
  ) + 
  geom_line(alpha = 0.5) +
  ggtitle("Zoom-in to GLM peaks") +
  theme_bw() +
  facet_wrap(variable~zoom_win, ncol = 2, scales = "free") ->
  pi_d_plot_zoom

ggsave(pi_d_plot_zoom, file = "pi_d_plot_zoom.pdf", h= 6, w = 8)






