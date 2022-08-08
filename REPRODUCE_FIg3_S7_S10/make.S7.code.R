### Make S7:
### 

library(vroom)
library(tidyverse)
library(magrittr)
library(reshape2)
library(foreach)

final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.12, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )

D.inv = vroom("Pi.TajD.genomewide/D.CM.W_100000.S_50000.INV.Tajima.D")
D.inv %<>% mutate(karyo = "inv")
D.std = vroom("Pi.TajD.genomewide/D.CM.W_100000.S_50000.STD.Tajima.D")
D.std %<>% mutate(karyo = "std")

rbind(D.inv, D.std) %>%
  as.data.frame() %>%
  filter(CHROM != "X") -> D.dat

D.dat %>% 
  filter(CHROM == "2L") %>%
  mutate(win =  case_when( BIN_START >= 6000000 & BIN_START <=  6400000 ~ "w6.1",
                            BIN_START >= 6600000 & BIN_START <=  7000000 ~ "w6.8",
                            BIN_START >= 9400000 & BIN_START <=  9800000 ~ "w9.6",
                            TRUE ~ "All.2L")) %>%
  group_by(win, karyo) %>%
  summarise(me = median(TajimaD)) -> D.values.win

## plot
  ggplot() +
  geom_density(data = D.dat, aes(
    x=TajimaD,
    color=karyo)) +
  geom_vline(data = D.values.win, aes(xintercept=me, color = karyo)) +
  geom_text(data = D.values.win, aes(x=me, y = 0,   label = win)) +
  facet_grid(~karyo) +
    theme_bw() +
    scale_color_manual(values = c("firebrick", "steelblue")) +
  geom_vline(xintercept = 0, linetype = "dashed") ->
    d.wins
  
  ggsave(d.wins, file = "d.wins.pdf", h = 3, w = 8)

#### some stats
#### Rankes value of 6.1
 sum(sort(D.dat$TajimaD[complete.cases(D.dat$TajimaD)]) < 
       filter(D.values.win, win == "w6.1", karyo == "inv" )$me) -> rank
 length(sort(D.dat$TajimaD[complete.cases(D.dat$TajimaD)])) -> all
 rank/all

 #### Rankes value of 9.6
 sum(sort(D.dat$TajimaD[complete.cases(D.dat$TajimaD)]) < 
       filter(D.values.win, win == "w9.6", karyo == "inv" )$me) -> rank
 length(sort(D.dat$TajimaD[complete.cases(D.dat$TajimaD)])) -> all
 rank/all
 
 ## general stats
 D.dat %>%
   group_by(karyo) %>%
   summarise(me = mean(TajimaD, na.rm = T),
             sd = sd(TajimaD, na.rm = T))
 
 
 ##### tmrca
 load("GEVA.2L.win.age.Rdata")
 
 final.windows.pos = 
   data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
              mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
   ) %>%
   mutate(start = (mid-0.2)*1e6 ,
          end  = (mid+0.2)*1e6  )
 
 ggplot() +
   geom_rect(data = final.windows.pos,
             aes(xmin=start, xmax = end,
                 ymin = 0, ymax = 190000), 
             alpha = 0.7, fill = "gold") +
   geom_line(
     data=filter(win.out.geva, V11 != "Dmel_ALL"),
     aes(
       x=mean.pos,
       y=(med.age*0.06666667),
       color = V11
     )) +
   #geom_ribbon() +
   geom_line() +
   ylab("TMRCA gens (x1e6)") +
   xlab("Genomic Position (Mb)") +
   ggtitle("TMRCA years") +
   theme_bw() +
   xlim(0,20.5*1e6) +
   geom_vline(xintercept = 2225744) +
   geom_vline(xintercept = 13154180) ->
   geva.win
 
 ggsave(geva.win, file ="geva.win.pdf", w = 7, h = 3)
 
 
 