library(tidyverse)
library(data.table)
library(magrittr)

final.windows.pos = 
  data.frame(win.name = c( "win_3.1", "win_4.7", "win_5.2", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )
setDT(final.windows.pos)
setkey(final.windows.pos, start, end)

win.age <- get(load("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/CODE/9.0.GEVA_allele_age/data/GEVA.2L.win.age.Rdata"))
win.age %>%
  filter(V11 == "Dmel_ALL") ->
  win.age.all

setDT(win.age.all)
win.age.all %<>%
  mutate(start = mean.pos, end =mean.pos)
setkey(win.age.all, start, end)


foverlaps(win.age.all,final.windows.pos) %>%
  mutate(Inversion = case_when(mean.pos < 2.2e6 & mean.pos > 13.2 ~ "outside",
                         TRUE ~ "inside")) %>%
  ggplot(aes(y=med.age*0.06666667, x =win.name, fill = Inversion)) +
  geom_boxplot() +
  geom_hline(yintercept = 85000) +
  coord_flip() +
  ylab("TRMCA (years)") + xlab("Window")

