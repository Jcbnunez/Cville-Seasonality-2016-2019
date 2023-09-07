### Visualize Analysis of New models

rm(list = ls())

### libraries
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

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))
sets

#### load models -->
t015 <- get(load("Window_analysis_temp.max;2;5.Cville.Rdata"))
t015 %>%
  mutate(data_type = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm")) %>%
  group_by(data_type, pos_mean, chr) %>%
  summarise(uci = quantile(rnp.binom.p, 0.01)) %>%
  mutate(model = "Tmax0-15") ->
  t015.ag

t015 %>%
  mutate(data_type = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm")) %>%
  group_by(data_type, pos_mean, chr) %>%
  summarise(uci.wza = quantile(rnp.wZa.p, 0.01)) %>%
  mutate(model = "Tmax0-15") ->
  t015.ag.wza




t715 <- get(load("Window_analysis_temp.max;4;5.Cville.Rdata"))
t715 %>%
  mutate(data_type = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm")) %>%
  group_by(data_type, pos_mean, chr) %>%
  summarise(uci = quantile(rnp.binom.p, 0.01)) %>%
  mutate(model = "Tmax7-15") ->
  t715.ag


EUW <- get(load("Window_analysis_humidity.ave;10;1.Europe_W.Rdata"))
EUW %>%
  mutate(data_type = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm")) %>%
  group_by(data_type, pos_mean, chr) %>%
  summarise(uci = quantile(rnp.binom.p, 0.01)) %>%
  mutate(model = "EUW") ->
  EUW.ag

 
EUE <- get(load("Window_analysis_humidity.var;3;3.Europe_E.Rdata"))
EUE %>%
  mutate(data_type = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm")) %>%
  group_by(data_type, pos_mean, chr) %>%
  summarise(uci = quantile(rnp.binom.p, 0.01)) %>%
  mutate(model = "EUE") ->
  EUE.ag

NAM <- get(load("Window_analysis_temp.min;8;2.North_America_E.Rdata"))
NAM %>%
  mutate(data_type = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm")) %>%
  group_by(data_type, pos_mean, chr) %>%
  summarise(uci = quantile(rnp.binom.p, 0.01)) %>%
  mutate(model = "NAM") ->
  NAM.ag

### PLOT Single models
#### Cville
 t015.ag %>%
 ggplot(aes(
 x=pos_mean/1e6,
 y=-log10(uci),
 color = data_type
 )) +
 geom_line() +
 facet_grid(.~chr, scale = "free_x") +
 theme_bw() ->
 t015.p.plot
ggsave(t015.p.plot, file = "t015.p.plot.pdf", h = 3, w = 9)

###
 rbind(
 t015.ag,
 t715.ag,
 EUW.ag,
 EUE.ag,
 NAM.ag
  ) %>%
 mutate(win_id = paste(chr, pos_mean, sep = "_")) ->
 tmods.ag

 rbind(
 EUW.ag,
 EUE.ag,
  NAM.ag ) %>%
 mutate(win_id = paste(chr, pos_mean, sep = "_")) ->
 Europa.ag


###
 tmods.ag %>%
 ggplot(aes(
 x=pos_mean/1e6,
 y=-log10(uci),
 color = data_type
 )) +
 geom_line() +
 facet_grid(model~chr, scale = "free_x") +
 theme_bw() ->
 all.p.plot
ggsave(all.p.plot, file = "all.p.plot.pdf", h = 8, w = 9)

####
 rbind(
tmods.ag ) %>%
 filter(data_type == "real") %>%
dcast(pos_mean+chr+win_id~model, value.var = "uci") ->
dcastdat

dcastdat %>%
dcast(pos_mean+chr~model, value.var = "uci") %>% 
ggplot(aes(
x=-log10(`Tmax0-15`),
y=-log10(`Tmax7-15`)
)) +
geom_point() +
theme_bw() +
 facet_wrap(.~chr, scale = "free", nrow = 1) +
geom_smooth(method = "lm") ->
 mods.cor
ggsave(mods.cor, file = "mods.cor.pdf", h = 3, w = 9)

###
### --> correlations
cor.test(~`Tmax0-15`+`Tmax7-15`, data = dcastdat)
cor.test(~`Tmax0-15`+`Tmax7-15`, data = filter(dcastdat, chr == "2L"))
cor.test(~`Tmax0-15`+`Tmax7-15`, data = filter(dcastdat, chr == "2R"))
cor.test(~`Tmax0-15`+`Tmax7-15`, data = filter(dcastdat, chr == "3L"))
cor.test(~`Tmax0-15`+`Tmax7-15`, data = filter(dcastdat, chr == "3R"))

###
dcastdat %>%
filter(chr == "2L") %>%
filter(`Tmax0-15` < 1e-60 & `Tmax7-15` < 1e-60) %>%
.$win_id -> pass.both

###
final.windows.pos = 
  data.frame(win.name = c(#"left", 
                          "win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6"
                          #, "right" 
                          ),
             mid = c(#2.2, 
                     3.1, 4.7, 5.2, 6.1, 6.8 , 9.6
                     #, 13.1
                     ),
             chr = "2L"
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )


### all pops
 tmods.ag %>%
 filter(chr == "2L") ->
 tmods.ag.2L
 
##  t015.ag.wza <------ i am here! 

 
 ggplot() +
   geom_rect(data = final.windows.pos, aes(xmin=start/1e6 , xmax =end/1e6, ymin = 0, ymax = 110),
              fill = "gold", alpha = 0.6) +
  geom_vline(xintercept = 2225744/1e6) +
  geom_vline(xintercept = 13154180/1e6) +
 geom_line(
 data = tmods.ag.2L,
 aes(
 x=pos_mean/1e6,
 y=-log10(uci),
 color = data_type
 )) +
 facet_grid(model~chr, scale = "free") +
 theme_bw() ->
 chr2L.p.plot
ggsave(chr2L.p.plot, file = "chr2L.p.plot.pdf", h = 6, w = 7)

### just charlottesville
### just charlottesville
### just charlottesville
### just charlottesville
### just charlottesville

###### ---> Must load t015 and process the two AG routines!

t015.ag %>%
filter(chr == "2L") ->
t015.ag.2L

left_join(
t015.ag.2L,
filter(t015.ag.wza, chr == "2L")
) -> t015.ag.2L.wza

 ggplot() +
   geom_rect(data = final.windows.pos, 
   aes(xmin=start/1e6 , xmax =end/1e6, 
   ymin = -110, ymax = 110),
              fill = "gold", alpha = 0.6) +
  geom_vline(xintercept = 2225744/1e6) +
  geom_vline(xintercept = 13154180/1e6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
 geom_line(
 data = t015.ag.2L.wza,
 aes(
 x=pos_mean/1e6,
 y=-log10(uci),
 color = data_type,
 linewidth = data_type
 )) +
  geom_line(
 data = t015.ag.2L.wza,
 aes(
 x=pos_mean/1e6,
 y=log10(uci.wza),
 color = data_type,
 linewidth = data_type
 )) +
scale_discrete_manual("linewidth", values = c(0.9, 1.2)) +
 xlim(0, 20.5) +
 facet_grid(model~chr, scale = "free") +
 theme_bw() ->
 chr2L.p.plot.cville
 
ggsave(chr2L.p.plot.cville, 
file = "chr2L.p.plot.cville.pdf", h = 2.3, w = 7)

###
###
###
###
###
###

###
  Europa.ag %>%
 filter(chr == "2L") ->
 tmods.ag.2L
 
 ggplot() +
   geom_rect(data = final.windows.pos, aes(xmin=start/1e6 , xmax =end/1e6, ymin = 0, ymax = 40),
              fill = "gold", alpha = 0.6) +
  geom_vline(xintercept = 2225744/1e6) +
  geom_vline(xintercept = 13154180/1e6) +
 geom_line(
 data = tmods.ag.2L,
 aes(
 x=pos_mean/1e6,
 y=-log10(uci),
 color = data_type
 )) +
 facet_grid(model~chr, scale = "free") +
 theme_bw() ->
 EU.chr2L.p.plot
ggsave(EU.chr2L.p.plot, file = "EU.chr2L.p.plot.pdf", h = 4, w = 7)

