## Explore LD correlations
## 
## 

library(data.table)
library(tidyverse)

dat <- fread("inv2L_informative_markers_Dm3.txt")


dat %>%
  .[order(as.numeric(pos)),] %>% 
  mutate(UCSC_format = paste("chr",chr,":",pos,"-",pos, sep = "") ) ->
  liftOver_positioninDGRP


write.table(liftOver_positioninDGRP$UCSC_format, 
            file = "./in2lt_prelift.txt", 
            append = FALSE, 
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = FALSE)

