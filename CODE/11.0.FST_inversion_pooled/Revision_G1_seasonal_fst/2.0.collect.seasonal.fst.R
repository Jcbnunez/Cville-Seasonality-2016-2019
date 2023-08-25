rm(list = ls())

library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(magrittr)
library(SeqArray)
library(lubridate)
library(gtools)
library(poolfstat)
library(geosphere)

final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )


root = "/gpfs2/scratch/jcnunez/genetics_resub/3.seasonal_fst/out_file/"

files.to.load = system(paste("ls ", root, "*", sep = ""), intern = T)

o <- 
foreach(fil = files.to.load,
.combine = "rbind")%do%{
tmp <- get(load(fil))
return(tmp)
}

###
o %>%
mutate(delta_T_group2 =
case_when(
Tdelta < 2 ~ "dT<2",
Tdelta > 10 ~ "dT>10",
)) %>% 
filter(delta_T_group2 %in% c("dT<2", "dT>10")) %>%
group_by(delta_T_group2, chr, start, end #, year
) %>%
summarize(
m.FST = mean(FST),
) -> dat.for.fig

ggplot() +
   geom_vline(xintercept = 2225744/1e6) +
   geom_vline(xintercept = 13154180/1e6) +
geom_rect(data = final.windows.pos,
            aes(xmin=start/1e6, xmax = end/1e6,
                ymin = -0.001, ymax = 0.014), 
            alpha = 0.7, fill = "gold") +
geom_line(data =dat.for.fig , aes(
x=(start+end)/(2*1e6),
y=m.FST,
color = delta_T_group2
),
linewidth = 1.05
) +
theme_bw() +
#facet_grid(year~.) + 
  xlim(0,20.5) +
  ylim(-0.001,0.014) ->
fst.plot.seas

ggsave(fst.plot.seas, 
file = "delta_T_group.pdf", 
w = 7, h = 2.4)

#### Some deltas
o %>%
mutate(delta_T_group2 =
case_when(
Tdelta < 2 ~ "dT<2",
Tdelta > 10 ~ "dT>10",
)) %>% 
  mutate(WinOverall = case_when(
    start > 2000000 & end < 2400000 ~ "bkpt",
    start > 12900000 & end < 13300000 ~ "bktp",
    start > 2800000 & end < 3200000 ~ "W",
    start > 4470000 & end < 4870000 ~ "W",
    start > 4920000 & end < 5320000 ~ "W",
    start > 6000000 & end < 6400000 ~ "W",
    start > 6600000 & end < 7000000 ~ "W",
    start > 9400000 & end < 9800000 ~ "W",
    TRUE ~ "NW")) %>% 
  mutate(WinSpec = case_when(
    start > 2000000 & end < 2400000 ~ "L",
    start > 12900000 & end < 13300000 ~ "R",
    start > 2800000 & end < 3200000 ~ "W3,1",
    start > 4470000 & end < 4870000 ~ "W4.6",
    start > 4920000 & end < 5320000 ~ "W5.2",
    start > 6000000 & end < 6400000 ~ "W6.1",
    start > 6600000 & end < 7000000 ~ "W6.8",
    start > 9400000 & end < 9800000 ~ "W9.6",
    TRUE ~ "NW")) ->
    o.annot

o.annot %>%
group_by(WinOverall, delta_T_group2#, WinSpec
) %>%
filter(!is.na(delta_T_group2)) %>%
summarize(m.FST = mean(FST)) 
 
   
summary(lm(FST~Tdelta, data = filter(o.annot, WinOverall == "NW" & delta_T_group2 == "dT>10")))

summary(lm(FST~Tdelta, data = filter(o.annot, WinOverall == "W" & delta_T_group2 == "dT>10")))

summary(lm(FST~Tdelta, data = filter(o.annot, WinSpec == "W5.2" & delta_T_group2 == "dT>10")))

o.annot %>%
filter(delta_T_group2 == "dT>10") %>%
group_by(WinSpec) %>%
slice_max(FST) %>% as.data.frame

#### PLot temperatures
o %>%
filter(!is.na(Tdelta) & !is.na(year)) %>%
group_by(samp1, year) %>%
summarize(m.T = mean(Temp1),
		m.dt = mean(Tdelta)) -> temps
		
temps %>%
ggplot(aes(
y=m.dt,
x=as.factor(year),
color = year
)) +
geom_boxplot() +
geom_point() ->
temps.raw

ggsave(temps.raw, file = "temps.raw.pdf")




#### Not used, yet interesting
o %>%
  mutate(Win = case_when(
    start > 2000000 & end < 2400000 ~ "left",
    start > 12900000 & end < 13300000 ~ "right",
    start > 2800000 & end < 3200000 ~ "w3.1",
    start > 4470000 & end < 4870000 ~ "w4.7",
    start > 4920000 & end < 5320000 ~ "w5.1",
    start > 6000000 & end < 6400000 ~ "w6.1",
    start > 6600000 & end < 7000000 ~ "w6.8",
    start > 9400000 & end < 9800000 ~ "w9.6",
    TRUE ~ "Other")) %>%
filter(!is.na(Tdelta) & !is.na(year)) %>%
filter(Win %in% c("Other", "w5.1") ) %>%
ggplot(aes(
x=as.factor(round(Tdelta, 0)),
y=FST,
)) +
geom_boxplot() +
facet_grid(Win~year) ->
fst.temp

ggsave(fst.temp, file = "fst.temp.pdf", w = 12, h = 5)



