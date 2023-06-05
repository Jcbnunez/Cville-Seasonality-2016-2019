## analysis of models

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

###
rm(list = ls())

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))
sets


###
results.folder <- "/scratch/yey2sn/old_scra/Overwintering_ms/4.2.env.omibus.mods/GLM_omnibus_window_analysis"

files <- system( paste("ls", results.folder, "| grep 'Win'  "), intern = T)

win.all.in =
  foreach(i=1:length(files), 
          .combine = "rbind")%do%{
            print(i)
            load(paste(results.folder, files[i], sep = "/"))
            win.out
          }

### Plot distributions

###
###
###EU-E (Tave45-75d; 9.4 times higher than permutations; Fig. S11A),
###
### average humidity 15-45 days prior for EU-W (H%ave15-45d; 5.3 times higher than permutations; Fig. S1
### 

final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )



win.all.in %>%
  filter(chr == "2L") %>%
  filter(model.pop == "temp.ave;9;3.Europe_E") %>%
  mutate(perm_type = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm",
  )) %>%
  group_by(chr, perm_type, pos_mean) %>%
  summarize(uci = quantile(rnp.binom.p, 0.05))  %>%
  mutate(metric = "rnvp") -> rnvpEUE
  
rnvpEUE %>%
  ggplot() +
  geom_rect(data = final.windows.pos,
            aes(xmin=start/1e6, xmax = end/1e6,
                ymin = -110, ymax = 110), 
            alpha = 0.7, fill = "gold") +
  geom_line(aes(
    x=pos_mean/1e6,
    y=-log10(uci),
    color=perm_type
  )) +
  geom_vline(xintercept = 2225744/1e6) +
  geom_vline(xintercept = 13154180/1e6) ->
  EUErnpv

ggsave(EUErnpv, file = "EUErnpv.pdf", h = 3, w = 9)

####
####

win.all.in %>%
  filter(chr == "2L") %>%
  filter(model.pop == "humidity.ave;8;1.Europe_W") %>%
  mutate(perm_type = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm",
  )) %>%
  group_by(chr, perm_type, pos_mean) %>%
  summarize(uci = quantile(rnp.binom.p, 0.05))  %>%
  mutate(metric = "rnvp") -> rnvpEUW

rnvpEUW %>%
  ggplot(aes(
    x=pos_mean/1e6,
    y=-log10(uci),
    color=perm_type
  )) +
  geom_line() +
  geom_vline(xintercept = 2225744/1e6) +
  geom_vline(xintercept = 13154180/1e6) ->
  EUWrnpv

ggsave(EUWrnpv, file = "EUWrnpv.pdf", h = 3, w = 9)

