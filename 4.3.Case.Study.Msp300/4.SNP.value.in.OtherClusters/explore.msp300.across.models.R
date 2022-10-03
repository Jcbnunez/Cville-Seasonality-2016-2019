#### How does the SNP perform in other contexts
#### 

library(tidyverse)
library(magrittr)
library(reshape2)

glm.cville <- get(load("/project/berglandlab/alan/environmental_ombibus_global/temp.max;2;5.Cville/temp.max;2;5.Cville.glmRNP.Rdata"))
glm.cville.0 = glm.cville %>%
  filter(perm == 0)

glm.EUE <- get(load("/project/berglandlab/alan/environmental_ombibus_global/temp.ave;9;3.Europe_E/temp.ave;9;3.Europe_E.glmRNP.Rdata"))
glm.EUE.0 = glm.EUE %>%
  filter(perm == 0)

glm.EUW <- get(load("/project/berglandlab/alan/environmental_ombibus_global/humidity.ave;8;1.Europe_W/humidity.ave;8;1.Europe_W.glmRNP.Rdata"))
glm.EUW.0 = glm.EUW %>%
  filter(perm == 0)

#"temp.ave;9;3.Europe_E",
#"humidity.ave;4;2.North_America_W",
#"humidity.ave;8;1.Europe_W"

glm.cville.0 %>%
  filter(chr == "2L" & pos == 5192177)
glm.EUE.0 %>%
  filter(chr == "2L" & pos == 5192177)
glm.EUW.0 %>%
  filter(chr == "2L" & pos == 5192177)





