library(tidyverse)
library(data.table)
args <- commandArgs(trailingOnly = T)
inprefix <- as.character(args[2])
bayesinfo <- read.table("/scratch/nzx3cc/nzx3cc/env_analysis/bayeswpos_1.out", header=T)[,1:2,]
chrlist <- c("2L", "2R", "3L", "3R")
preparer <- function(i){
  file <- read.table(paste0(i, "_", inprefix, "_summary_betai_reg.out"), header=T)
  relevant <- bayesinfo %>% filter(CHR == i)
  relevant <- rep(relevant, 10)
  fixedtab <- data.table(relevant, file)
  assign(paste0(i, "_table"), fixedtab)
}
tabs <- rbindlist(lapply(lapply(chrlist, preparer), as.data.frame.list))
tabs <- sort(tabs, CHR, POS)
write.table(tabs, paste0("LOCO_", inprefix, "bayes.txt"))
