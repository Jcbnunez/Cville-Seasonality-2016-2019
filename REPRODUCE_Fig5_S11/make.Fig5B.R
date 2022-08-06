### Make Panel 5B
### 

library(tidyverse)
library(magrittr)
library(data.table)
library(car)
library(DescTools)
library(foreach)
library(doMC)
library(patchwork)
library(ggbeeswarm)
library(reshape2)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(viridis)
library(lubridate)
library(forcats)
library(viridis)
library(SeqArray)
library(tidyverse)
library(gmodels)
library(scatterpie)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(weathermetrics)
registerDoMC(2)
library(broom)

####

load("DatFor.Haplotypes.trajectory.time.weather.Rdata")

###


#####
Cville_haplotags_for_viz %>%
  ggplot(aes(
    x=yday,
    y=Mean_haplotag,
    #ymin=ci_l,
    #ymax=ci_h,
    color=(temp.max),
  )) + 
  #geom_smooth(method = "lm", se = F, size = 0.8, color = "grey") +
  #geom_errorbar(width = 0.1) +
  scale_color_gradient2(low="steelblue", high = "firebrick2", mid = "gold1", 
                        midpoint = 25) +
  geom_point(aes(shape=as.factor(year))) +
  geom_smooth( se = F, size = 0.8, color = "black", linetype = "dashed") +
  ylim(0,0.38) +
  theme_bw() + 
  facet_grid(win~.)->
  haplo.time.colortemp.ave

ggsave(haplo.time.colortemp.ave, file ="haplo.time.colortemp.ave.pdf", h = 6, w = 3.2)


###
###


Cville_haplotags_for_viz %>%
  melt(id = c("sampleId", 
              "collectionDate", 
              "set", 
              "year", 
              "win", 
              "yday", 
              "Mean_haplotag")) %>% 
  filter(variable == "temp.max") %>%
  ggplot(aes(
    x=value,
    y=Mean_haplotag
  )) +
  geom_point(color = "grey",aes(shape=as.factor(year))) +
  geom_smooth(method = "lm", 
              color = "black") +
  theme_bw() +
  facet_grid(win~., scales = "free_x") ->
  eco.vars.afs

ggsave(eco.vars.afs, file ="eco.vars.afs.pdf", h = 6, w = 3.2)


### P-values for regressions:
### 

Cville_haplotags_for_viz %>%
  melt(id = c("sampleId", 
              "collectionDate", 
              "set", 
              "year", 
              "win", 
              "yday", 
              "Mean_haplotag")) %>% 
  filter(variable == "temp.max") %>% 
  nest(data = -c(win) ) %>% 
  mutate(model = map(data, ~lm(Mean_haplotag~value, data = .)),
         tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  filter(term == "value")  %>%
  dplyr::select(win, estimate, p.value, std.error) %>% 
  mutate(estimate = as.numeric(format(estimate, scientific = T, digits = 1)),
         p.value = as.numeric(format(p.value, scientific = T, digits = 1))) %>% 
  as.data.frame()


