### Examine new models for Genetics revision #1

library(tidyverse)
library(data.table)
library(magrittr)

#### object 

  load(newmods)
  o2.ag <- o2.ag[!cluster%in%c("2.North_America_Mid", "2.North_America_W")]
o2.ag <- o2.ag %>% filter(perm_strategy == "new")
  o2.ag[,sig:=log2(prop.real/prop.perm.uci)>0 | log2(prop.real/prop.perm.lci)<0]

#
###
