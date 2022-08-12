### Clear data

library(tidyverse)
library(reshape2)
library(magrittr)

Sechelia.uncleared <- read.delim("./Sechelia.uncleared.txt", na.strings = "")

Sechelia.uncleared %>% 
  .[complete.cases(.),] %>%
  group_by(samp) %>%
  slice_head() -> unq.schel

write.table(unq.schel$SRR, file = "Sechelia.SRRs.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(unq.schel, file = "Sechelia.metadata.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = T, qmethod = c("escape", "double"),
            fileEncoding = "")

#####
#####
#####

yakuba.uncleared <- read.delim("./yakuba.uncleared.txt", na.strings = "")

yakuba.uncleared %>% 
  .[complete.cases(.),] %>%
  group_by(samp) %>%
  slice_head() -> unq.yak

write.table(unq.yak$SRR, file = "yakuba.SRRs.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(unq.yak, file = "yakuba.metadata.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = T, qmethod = c("escape", "double"),
            fileEncoding = "")

####

write.table(c(unq.yak$SRR,unq.schel$SRR) , file = "yak.sec.SRAs.joint.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")


