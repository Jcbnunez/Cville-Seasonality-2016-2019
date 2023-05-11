### LD plot

library(tidyverse)
library(data.table)
library(magrittr)

setwd("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/CODE/8.0.Linkage/LD_of_Windows_May2023/")

###
ld.out <- get(load("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/CODE/8.0.Linkage/LD_of_Windows_May2023/ld_df_winanot_pair.flt.glm.scored.Rdata"))

setDT(ld.out)

ld.out %>%
  filter(!is.na(p_lrt.tresh)) %>%
  group_by(bkpt) %>%
  summarise(mean(R2))

ld.out %>%
  filter(!is.na(p_lrt.tresh)) %>%
  group_by(p_lrt.tresh) %>%
  summarise(mean(R2))


ld.out %>%
  filter(!is.na(p_lrt.tresh)) %>%
  ggplot(aes(
    x=outloc,
    y=R2,
    fill=p_lrt.tresh
  )) +
  facet_wrap(~bkpt) +
  theme_bw() +
  geom_boxplot(outlier.shape = NA) ->
  ld.outs

ggsave(ld.outs, file = "ld.outs.pdf", w = 7, h = 3)

####
load("GEVA.win.age.Rdata")
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


  ggplot() +
    geom_rect(data = final.windows.pos, 
              aes(xmin=start/1e6 , 
                  xmax =end/1e6, ymin = 0, 
                  ymax = 170000),
              fill = "gold", alpha = 0.6) +
  geom_line(linewidth = 0.7, color = "grey30",
  data=filter(win.out.geva, V11 %in% c("Dmel_het") ),
  aes(
    x=mean.pos/1e6,
    y=med.age*0.06666667,
    #color =V11
  )) +
    theme_bw() +
  geom_vline(xintercept = 2225744/1e6) +
  geom_vline(xintercept = 13154180/1e6) +
  geom_hline(yintercept = 95000, 
             linetype = "dashed", 
             color = "red")  ->
    age.inv
  
  ggsave(age.inv, file = "age.inv.pdf", w = 8, h = 2.3)
  
  
