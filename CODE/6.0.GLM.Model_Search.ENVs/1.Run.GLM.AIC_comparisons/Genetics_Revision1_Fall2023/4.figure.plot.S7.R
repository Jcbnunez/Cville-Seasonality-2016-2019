library(data.table)
library(tidyverse)

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

sets %>%
  ggplot(aes(
    x = -start,
    xend = -end,
    y = mod, 
    yend = mod
  )) +
  geom_segment(size = 1.2) +
  theme_classic()

all.models.info <- read.delim("~/Documents/GitHub/Cville-Seasonality-2016-2019/CODE/6.0.GLM.Model_Search.ENVs/1.Run.GLM.AIC_comparisons/Data/all.models.info.txt", header=FALSE)

all.models.info$V1 %>% unique()