### Prepare Fst sites for analysis
### 
### 

rm(list = ls())

library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(magrittr)
library(doMC)
registerDoMC(4)
library(SeqArray)
library(lubridate)
library(tidyverse)
library(vroom)

######

marker.metadat <- vroom("SNPs.for.fst.analysis.Metadat.txt")
names(marker.metadat)[2:3] = c("CHROM","POS")
fst.markers <- vroom("MARKERS.inv.vs.std.cm.fst.weir.fst")

left_join(fst.markers, marker.metadat) -> fst.metadat

fst.metadat$type[which(fst.metadat$type == "glm.mrks")] = "inv.markers"

fst.metadat %>%
  ggplot(aes(
    x=type,
    y=WEIR_AND_COCKERHAM_FST
  )) +
  geom_boxplot(outlier.shape =  NA) +
  theme_bw() +
  xlab("Marker type") +
  ylab(expression(F[ST]))->
  fst.box
ggsave(fst.box, file = "fst.box.pdf", w = 4, h = 3)
ggsave(fst.box, file = "fst.box.png", w = 4, h = 3)







