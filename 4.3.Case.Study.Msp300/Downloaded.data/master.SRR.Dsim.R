### Clear data

library(tidyverse)
library(reshape2)
library(magrittr)

Sechelia.uncleared <- read.delim("~/Documents/GitHub/Cville-Seasonality-2016-2019/4.3.Case.Study.Msp300/Downloaded.data/Sechelia.uncleared.txt")

data.frame(
samp= unique(Sechelia.uncleared$Sample),
SRR = unique(Sechelia.uncleared$SRR)) 
  