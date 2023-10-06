library(tidyverse)
library(magrittr)
library(vroom)

library(poolfstat)

library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(gtools)
library(patchwork)
library(foreach)
library(doMC)
library(doParallel)
registerDoMC(10)


args <- commandArgs(trailingOnly = TRUE)
file1 <- as.character(args[2])
cov <- as.numeric(args[3])
bayes <- read.table(file1, header=T)
wins <- read.table("windows.out", header=T)

bucket <- data.table(c(1))
> for (j in 1:1925){
+ tmp <- wins %>% filter(i == j)
+ bayeswind <- bayes %>% filter(POS > tmp$start, POS > tmp$end, COVARIABLE==cov)
+ bayeswindlist <- bayeswind %>% filter(eBPis > 1.3)
+ success <- (dim(bayeswindlist)[1])
+ total <- (dim(bayeswind)[1])
+ pval <- binom.test(success, total, 0.05)
+ bucket <- merge(bucket, pval, all=TRUE)
 }
bucket <- bucket[2:1925,]
write.table(bucket, "binomwindow.out", row.names=F, col.names=T)


