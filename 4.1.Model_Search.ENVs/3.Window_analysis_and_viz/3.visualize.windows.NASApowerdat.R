#### Graph models
#### 
## analysis of models
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

results.folder <- "/scratch/yey2sn/Overwintering_ms/4.2.env.omibus.mods/GLM_omnibus_window_analysis/0.01"

files <- system( paste("ls", results.folder, "| grep 'Win'  "), intern = T)

win.all.in =
foreach(i=1:length(files), .combine = "rbind")%do%{
  print(i)
  load(paste(results.folder, files[i], sep = "/"))
  win.out
}

##----> visualize
data.frame(
  start= c(4650065, 5100324, 6100321, 9500286),
  end= c(4799922, 5349218, 6349489, 9700005)
) %>% mutate(mid.point = (start/2) + (end/2) ) %>%
  mutate(win.name = paste("win", round(mid.point/1e6, 1), sep = "."  )) ->
  final.windows.pos

win.all.in %>%
   mutate(data_type = case_when(perm == 0 ~ "real",
                                perm != 0 ~ "perm")) %>%
   group_by(data_type, pos_mean, model.pop) %>%
   summarise(uci = quantile(rnp.binom.p, 0.05)) ->
    dat.for.plot
  
dat.for.plot %>% 
  dcast(pos_mean+model.pop~data_type) %>% 
  filter(-log10(real) > -log10(perm)) %>%
  mutate(uci = real) ->
  dat.for.plot.sig

dat.for.plot.sig%>%
  dplyr::select(pos_mean,model.pop, uci) %>%
  dcast(pos_mean~model.pop)  %>% 
  mutate(na.count = rowSums(is.na(.))) %>% 
  filter(na.count <= 3) %>%
  .$pos -> sig.pos

#####
   ggplot() +
   geom_line(
     data=dat.for.plot,
     aes(
       x=pos_mean/1e6,
       y=-log10(uci),
       color = data_type
     ), alpha = 0.5
   ) +
  geom_rect(data = final.windows.pos,
            aes(xmin=start/1e6, xmax = end/1e6,
                ymin = -0, ymax = 45), 
            alpha = 0.7, fill = "gold") +
     geom_point(
       data=filter(dat.for.plot.sig, pos_mean %in% sig.pos),
       aes(
         x=pos_mean/1e6,
         y=-log10(uci),
        
         #color = data_type
       ),  size = 1.4
     ) +
   geom_vline(xintercept = 2225744/1e6) +
   geom_vline(xintercept = 13154180/1e6) +
   facet_grid(model.pop~.)->
   tester_fig
 
 ggsave(tester_fig, file = "tester_fig.pdf", w = 12, h = 8)

 
 
 ######
 ######
 ######
 

 