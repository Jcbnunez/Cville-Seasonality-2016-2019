#inour et al 1984
library(tidyverse)
library(magrittr)

inoue.1984 <- read.delim("~/Documents/GitHub/Cville-Seasonality-2016-2019/20.old.papers.dat/inoue.1984.txt")

inoue.1984 %>% 
  ggplot(aes(
    x=Month,
    y=In2Lt,
    color=as.factor(Year)
  )) +
  geom_smooth(se =F) +
  geom_point()