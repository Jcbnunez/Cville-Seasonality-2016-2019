### Plot resampling
### 
rm(list = ls())

library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(gtools)
library(patchwork)
library(ggalt)
library(reshape2)
library(ggrepel)
library(sp)
library(permute)
library(RColorBrewer)
library(devtools)
library(lubridate)
#install_github('tavareshugo/windowscanr')
library(windowscanr)


load("./remaple_lists.Rdata")


real_df = do.call(rbind, loop2_list)
