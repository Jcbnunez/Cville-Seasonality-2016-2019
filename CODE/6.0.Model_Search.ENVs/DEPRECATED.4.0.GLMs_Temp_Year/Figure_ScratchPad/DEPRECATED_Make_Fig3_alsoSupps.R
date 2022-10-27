### Prepare Panels for figure 4. 
### 

library(tidyverse)
library(vroom)
library(tidyverse)
library(data.table)
library(tidyr)
library(viridis)
library(patchwork)
library(magrittr)
library(RColorBrewer)
library(rstatix)


#### Part 1 --- Describe the distribution of p.values

### ---> SEE special script ---> 

##### Temperature Figure -- Not part of figure 3

#github_addr <- "/Users/jcbnunez/Documents/GitHub"
### Part 1: load data
#load(paste(github_addr,"/Cville-Seasonality-2016-2019/4.GLMs_Temp_Year/temperatureAverage_yearFactor_GLM/weatherAve#.Rdata", sep =""))
#
#metadat <- vroom(paste(github_addr,"/Cville-Seasonality-2016-2019/4.GLMs_Temp_Year/temperatureAverage_yearFactor_GLM#/DEST_10Mar2021_POP_metadata.csv", sep =""))
#
#names(weather.ave)[1] = "sampleId"
#
#weather.ave %>%
#  left_join(metadat) %>%
#  mutate(date_posix = as.Date(collectionDate, format = "%m/%d/%Y")) %>%
#  ggplot(aes(x = date_posix, 
#             y = aveTemp/10,
#             group = year,
#             fill = aveTemp/10
#             )) +
#  geom_line() +
#  geom_point(shape = 21, size = 1.5) +
#  scale_x_date(date_labels = "%y") +
#  scale_fill_gradient2(low = "steelblue", high = "firebrick", midpoint = 15) +
#  theme_bw() +
#  theme(legend.position = "top") +
#  facet_grid(.~locality, scales = "free_x",  space = "free_x")
#
#########
######### Plot the giant manhattan plot
#

