library(tidyverse)
library(reshape2)
library(magrittr)

samples_to_phase <- read.delim("~/Dropbox/2021_Cville_DEST/Alyssa's Single Ind/Guide_Files/samples_to_phase.txt", header=FALSE)

samples_to_phase %<>%
  mutate(`2L` = "chr",
         `2R` = "chr",
         `3L` = "chr",
         `3R` = "chr",
         `X` = "chr",
         ) %>% 
  melt(id = c("V1","V2"))

write.table(samples_to_phase[,c("V1","V2","variable")],
            sep = "\t",
            row.names = F,
            col.names = F,
            quote = F,
            file = "~/Dropbox/2021_Cville_DEST/Alyssa's Single Ind/Guide_Files/samples_to_phase_chr.txt")

