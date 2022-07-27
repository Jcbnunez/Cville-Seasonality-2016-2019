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


load("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/nasa_power.weather.mod.Rdata")
names(weather.ave)[1] = "sampleId"

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

weather.ave %>%
  filter(mod == 2) %>%
  dplyr::select(sampleId, temp.max) -> tmax

weather.ave %>%
  filter(mod == 9) %>%
  dplyr::select(sampleId, temp.ave) -> tave

weather.ave %>%
  filter(mod == 8) %>%
  dplyr::select(sampleId, humidity.ave) -> have


cbind(tmax,  tave[,-1], have[,-1]) %>%
  as.data.frame() %>%
  group_by(sampleId) %>%
  slice_head(n=1) -> obs_for_cors

cor.test(~ temp.max + temp.ave, dat = obs_for_cors)
cor.test(~ temp.max + humidity.ave, dat = obs_for_cors)
cor.test(~ temp.ave + humidity.ave, dat = obs_for_cors)

